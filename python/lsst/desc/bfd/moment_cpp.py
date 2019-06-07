import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin
from lsst.meas.base.pluginRegistry import register
import numpy as np
from scipy.optimize import fsolve

import bfd
from .config import BFDConfig
from . import KSigmaWeightF
__all__ = (
    "BfdMomentConfig", "BfdMoment"
)


class CenterMoment():

    def __init__(self, bfd_config, kgal):
        self.bfd_config = bfd_config
        self.kgal = kgal

    def xy_moment(self, dx):
        '''Interface to fsolve to return x,y moments given input origin shift
            '''
        return self.kgal.getShifted(dx[0], dx[1]).getMoment().xy

    def xy_jacobian(self, dx):
        '''Function to return Jacobian of X & Y moments with respect to
        origin shift dx, for use in solver.  Use the fact that posn derivatives
        of the first moments are the second moments.
        '''
        mom = self.kgal.getShifted(dx[0], dx[1]).getMoment()
        e = mom.m
        return -0.5 * np.array([[e[self.bfd_config.MR] + e[self.bfd_config.M1],
                                 e[self.bfd_config.M2]],
                                [e[self.bfd_config.M2],
                                 e[self.bfd_config.MR] - e[self.bfd_config.M1]]])

    def recenter(self, sigma):
        '''Find dx, dy that null the X and Y moments

        Returns:
        dx    Shift of origin applied
        error True if there was a failure to converge
        msg   message string on failure
        '''
        dx = np.zeros(2, dtype=float)
        dx, junk, ier, msg = fsolve(self.xy_moment, dx, fprime=self.xy_jacobian, full_output=True)

        threshold = np.sqrt(2.0) * np.sqrt(sigma**2)
        wandered_too_far = np.abs(dx) >= threshold
        badcentering = wandered_too_far[0] or wandered_too_far[1] or ier <= 0
        return dx, badcentering, msg


class BfdMomentConfig(SingleFramePluginConfig):
    use_mag = pexConfig.Field(dtype=bool, default=True, doc="include magnification")
    use_conc = pexConfig.Field(dtype=bool, default=True, doc="include concentration")
    weight_n = pexConfig.Field(dtype=int, default=4, doc="n for k-sigma weight function")
    weight_sigma = pexConfig.Field(dtype=float, default=1.25, doc="sigma for k-sigma weight function")
    max_shift = pexConfig.Field(dtype=float, default=4, doc="maximum centroid shift")
    min_size = pexConfig.Field(dtype=int, default=48, doc="minimum image size")


@register("BfdMoment")
class BfdMoment(SingleFramePlugin):
    '''
    Compute the moments of an image using the k-sigma weight function.  
    '''

    ConfigClass = BfdMomentConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.bfd = BFDConfig(use_conc=self.config.use_conc, use_mag=self.config.use_mag)
        self.n_even = self.bfd.BFDConfig.MSIZE
        self.n_odd = self.bfd.BFDConfig.XYSIZE

        self.moment = schema.addField('bfd_moments', type="ArrayF",
                                      size=self.n_even, doc="Even Bfd moments")
        self.failed_moment = schema.addField('failed_moments', type="ArrayF",
                                      size=self.n_even, doc="Even Bfd moments")
        self.odd = schema.addField('bfd_odd', type="ArrayF",
                                   size=self.n_odd, doc="odd moments")
        self.failed_odd = schema.addField('failed_odd', type="ArrayF",
                                   size=self.n_odd, doc="odd moments")
        self.shift = schema.addField('bfd_shift', type="ArrayF",
                                     size=2, doc="amount shifted to null moments")
        self.cov_even = schema.addField('bfd_cov_even', type="ArrayF",
                                        size=self.n_even*(self.n_even+1)//2,
                                        doc="even moment covariance matrix")
        self.cov_odd = schema.addField('bfd_cov_odd', type="ArrayF",
                                       size=self.n_odd*(self.n_odd+1)//2,
                                       doc="odd moment covariance matrix")
        self.flag = schema.addField('bfd_flag', type="Flag", doc="Set to 1 for any fatal failure")
        self.centroid_flag = schema.addField('bfd_flag_centroid', type="Flag",
                                             doc="Set to 1 for any fatal failure of centroid")

        self.weight = KSigmaWeightF(self.config.weight_sigma, self.config.weight_n)

    def measure(self, measRecord, exposure):
        center = measRecord.getCentroid()
        local_lin_wcs = exposure.getWcs().linearizePixelToSky(center, afwGeom.arcseconds)

        jacobian = local_lin_wcs.getLinear().getMatrix()
        sky_pos = exposure.getWcs().pixelToSky(center)
        uvref = (sky_pos.getRa().asArcseconds(), sky_pos.getDec().asArcseconds())

        box = measRecord.getFootprint().getBBox()
        xy_pos = (center.getX() - box.getMinX(), center.getY() - box.getMinY())

        bfd_wcs = bfd.WCS(jacobian, xyref=xy_pos, uvref=uvref)

        factor = exposure.getMetadata().get('variance_scale')
        noise = np.median(exposure.variance[box].array)/factor

        image = exposure.image[box].array

        psf_image = exposure.getPsf().computeKernelImage(center).array

        kdata = bfd.generalImage(image, uvref, psf_image, wcs=bfd_wcs, pixel_noise=noise,
                                 size=self.config.min_size)
        conjugate = set(np.where(kdata.conjugate.flatten()==False)[0])
        kgal = self.bfd.KGalaxy(self.weight, kdata.kval.flatten(), kdata.kx.flatten(), kdata.ky.flatten(),
                                kdata.kvar.flatten(), kdata.d2k, conjugate)

        #nulled_data, dx, dy = kgal.getNulled(self.config.max_shift)

        cm = CenterMoment(self.bfd.BFDConfig, kgal)
        dx, failed, msg = cm.recenter(self.config.weight_sigma)
        #print(dx,converge,msg)
        if failed is True:
            target = kgal.getTarget()
            measRecord.set(self.flag, 1)
            measRecord.set(self.centroid_flag, 1)
            measRecord.set(self.failed_moments, np.array(target.mom.m, dtype=np.float32))
            measRecord.set(self.failed_odd, np.array(target.mom.xy, dtype=np.float32))
        else:
            target = kgal.getShifted(dx[0], dx[1]).getTarget()
            mom_even = target.mom.m
            mom_odd = target.mom.xy
            cov_even = target.cov.m
            cov_odd = target.cov.xy

            cov_even_save = []
            cov_odd_save = []
            for ii in range(cov_even.shape[0]):
                cov_even_save.extend(cov_even[ii][ii:])
            for ii in range(cov_odd.shape[0]):
                cov_odd_save.extend(cov_odd[ii][ii:])

            measRecord.set(self.moment, np.array(mom_even, dtype=np.float32))
            measRecord.set(self.odd, np.array(mom_odd, dtype=np.float32))
            measRecord.set(self.cov_even, np.array(cov_even_save, dtype=np.float32))
            measRecord.set(self.cov_odd, np.array(cov_odd_save, dtype=np.float32))
            measRecord.set(self.shift, np.array([dx[0], dx[1]], dtype=np.float32))

    def fail(self, measRecord, error=None):
        measRecord.set(self.flag, True)
