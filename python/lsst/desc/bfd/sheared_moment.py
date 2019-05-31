
import lsst.pex.config as pexConfig
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.geom as geom
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin
from lsst.meas.base.pluginRegistry import register
#from lsst.meas.base.baseMeasurement import BaseMeasurementPluginConfig
import numpy as np

import bfd
from .config import BFDConfig
import galsim

usingHSC = True
try:
    import hsc.pipe
except:
    usingHSC = False

__all__ = (
    "PythonBfdShearedMomentConfig", "PythonBfdShearedMoment"
    )
    

class PythonBfdShearedMomentConfig(SingleFramePluginConfig):
    weight_n = pexConfig.Field(dtype=int, default=4, doc="n for k-sigma weight function")
    weight_sigma = pexConfig.Field(dtype=float, default=1.25, doc="sigma for k-sigma weight function")
    max_g = pexConfig.Field(dtype=float, default=0.15, doc="maximum value of g")
    max_mu = pexConfig.Field(dtype=float, default=0.15, doc="min/max value of mu")
    sigma_dx = pexConfig.Field(dtype=float, default=1, doc="sigma of x,y in pixels")


@register("pythonBfdShearedMoment")
class PythonBfdShearedMoment(SingleFramePlugin):
    '''
    Compute the moments of an image using the k-sigma weight function.  This uses the new python BFD
    implementation.
    '''

    ConfigClass = PythonBfdShearedMomentConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)


        self.bfd_config = BFDConfig(use_mag=True, use_conc=True)

        NMoment = self.bfd_config.MSIZE
        self.moment = schema.addField('moment', type="ArrayF", size=NMoment, doc="Even Bfd moments")
        self.sheared_moment = schema.addField('sheared_moment', type="ArrayF", size=NMoment, doc="Sheared Even Bfd moments")
        self.flag = schema.addField('bfd.flags', type="Flag", doc="Set to 1 for any fatal failure")
        self.flag_sat = schema.addField('bfd.flags.saturated.center', type="Flag", doc="Set to 1 for any fatal failure")
        self.flag_parent = schema.addField('bfd.flags.parent', type="Flag", doc="Set to 1 for any fatal failure")
        self.centroid_flag = schema.addField('bfd_flag_centroid', type="Flag", doc="Set to 1 for any fatal failure of centroid")

        self.g1 = schema.addField('g1', type=float, doc="g1")
        self.g2 = schema.addField('g2', type=float, doc="g2")
        self.mu = schema.addField('mu', type=float, doc="mu")
        self.dx = schema.addField('dx', type=float, doc="dx")
        self.dy = schema.addField('dy', type=float, doc="dy")

        self.weight = bfd.KSigmaWeightF(float(self.config.weight_sigma), int(self.config.weight_n))

    def measure(self, measRecord, exposure):

        center = measRecord.getCentroid()
        psf_image = exposure.getPsf().computeKernelImage(center).array

        ra = measRecord.get('coord_ra').asArcseconds()
        dec = measRecord.get('coord_dec').asArcseconds()

        if usingHSC:
            local_lin_wcs = exposure.getWcs().linearizePixelToSky(center, afwGeom.arcseconds)
        else:
            local_lin_wcs = exposure.getWcs().linearizePixelToSky(center, geom.arcseconds)

        jacobian = local_lin_wcs.getLinear().getMatrix()
        sky_pos = exposure.getWcs().pixelToSky(center)
        uvref = (sky_pos.getRa().asArcseconds(), sky_pos.getDec().asArcseconds())

        if usingHSC:
            box = measRecord.getFootprint().getBBox()
            xy_pos = (center.getX() - box.getMinX(), center.getY() - box.getMinY())
            box.shift(-afwGeom.Extent2I(exposure.getXY0()))
        else:
            box = measRecord.getFootprint().getBBox(afwImage.PARENT)
            xy_pos = (center.getX() - box.getMinX(), center.getY() - box.getMinY())
        
        
        bfd_wcs = bfd.WCS(jacobian, xyref=xy_pos, uvref=uvref)
        image = exposure.image[box].array

        kdata = bfd.generalImage(image, uvref, psf_image, wcs=bfd_wcs, pixel_noise=1)
        conjugate = set(np.where(kdata.conjugate.flatten()==False)[0])

        kgal = self.bfd_config.KGalaxy(self.weight, kdata.kval.flatten(), kdata.kx.flatten(), kdata.ky.flatten(), 
                                      kdata.kvar.flatten(), kdata.d2k, conjugate);

        g1 = np.random.rand()*(2*self.config.max_g) - self.config.max_g
        g2 = np.random.rand()*(2*self.config.max_g) - self.config.max_g
        mu = np.random.rand()*(2*self.config.max_mu) - self.config.max_mu
        shift_dx = np.random.randn()*self.config.sigma_dx
        shift_dy = np.random.randn()*self.config.sigma_dx

        nulled_data, dx, dy = kgal.getNulled(4)
        if nulled_data is None:
            measRecord.set(self.flag, True)
            measRecord.set(self.centroid_flag, True)
            return
        
        target = nulled_data.getTarget()
        shifted_data = nulled_data.getShifted(shift_dx, shift_dy)
        sheared_target = shifted_data.getShearedTarget(g1, g2, mu, False)

        moment = target.mom.m
        sheared_moment = sheared_target.mom.m

        measRecord.set(self.moment, np.array(moment, dtype=np.float32))
        measRecord.set(self.sheared_moment, np.array(sheared_moment, dtype=np.float32))
        measRecord.set(self.g1, g1)
        measRecord.set(self.g2, g2)
        measRecord.set(self.mu, mu)
        measRecord.set(self.dx, shift_dx)
        measRecord.set(self.dy, shift_dy)

    def fail(self, measRecord, error=None):
        print('failed %s',error)
        measRecord.set(self.flag, True)

