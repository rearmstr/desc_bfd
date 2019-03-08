import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin
from lsst.meas.base.pluginRegistry import register
import numpy as np

import bfd

__all__ = (
    "PythonBfdMomentConfig", "PythonBfdMoment"
)

NMoment = 5


class PythonBfdMomentConfig(SingleFramePluginConfig):
    n_moment = pexConfig.Field(dtype=float, default=NMoment, doc="number of moments measured")
    weight_n = pexConfig.Field(dtype=float, default=4, doc="n for k-sigma weight function")
    weight_sigma = pexConfig.Field(dtype=float, default=1.25, doc="sigma for k-sigma weight function")


@register("pythonBfdMoment")
class PythonBfdMoment(SingleFramePlugin):
    '''
    Compute the moments of an image using the k-sigma weight function.  This uses the new python BFD
    implementation.
    '''

    ConfigClass = PythonBfdMomentConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.moment = schema.addField('bfd_moments', type="ArrayF", size=NMoment, doc="Even Bfd moments")
        self.xy = schema.addField('bfd_xy', type="ArrayF", size=2, doc="xy shift")
        self.cov_even = schema.addField('bfd_cov_even', type="ArrayF", size=NMoment*(NMoment+1)//2,
                                        doc="even moment covariance matrix")
        self.cov_odd = schema.addField('bfd_cov_odd', type="ArrayF", size=3,
                                       doc="odd moment covariance matrix")
        self.flag = schema.addField('bfd_flag', type="Flag", doc="Set to 1 for any fatal failure")
        self.centroid_flag = schema.addField('bfd_flag_centroid', type="Flag",
                                             doc="Set to 1 for any fatal failure of centroid")

        self.weight = bfd.KSigmaWeight(self.config.weight_n, self.config.weight_sigma)

    def measure(self, measRecord, exposure):
        center = measRecord.getCentroid()

        local_lin_wcs = exposure.getWcs().linearizePixelToSky(center, afwGeom.arcseconds)

        jacobian = local_lin_wcs.getLinear().getMatrix()
        sky_pos = exposure.getWcs().pixelToSky(center)
        uvref = (sky_pos.getRa().asArcseconds(), sky_pos.getDec().asArcseconds())

        box = measRecord.getFootprint().getBBox()
        xy_pos = (center.getX() - box.getMinX(), center.getY() - box.getMinY())

        box.shift(-afwGeom.Extent2I(exposure.getXY0()))
        bfd_wcs = bfd.WCS(jacobian, xyref=xy_pos, uvref=uvref)
        noise = np.median(exposure.variance[box].array)

        image = exposure.image[box].array

        psf_image = exposure.getPsf().computeKernelImage(center).array

        kdata = bfd.generalImage(image, uvref, psf_image, wcs=bfd_wcs, pixel_noise=noise)
        moment_calc = bfd.MomentCalculator(kdata, self.weight, id=measRecord.getId())

        xyshift, error, msg = moment_calc.recenter()
        if error:
            measRecord.set(self.flag, 1)
            measRecord.set(self.centroid_flag, 1)
            return
        else:
            covgal = moment_calc.get_covariance()
            moment = moment_calc.get_moment(0, 0)
            measRecord.set(self.moment, np.array(moment.even, dtype=np.float32))

            cov_even_save = []
            cov_odd_save = []
            for ii in range(covgal[0].shape[0]):
                cov_even_save.extend(covgal[0][ii][ii:])
            for ii in range(covgal[1].shape[0]):
                cov_odd_save.extend(covgal[1][ii][ii:])
            measRecord.set(self.cov_even, np.array(cov_even_save, dtype=np.float32))
            measRecord.set(self.cov_odd, np.array(cov_odd_save, dtype=np.float32))
            measRecord.set(self.xy, np.array(xyshift, dtype=np.float32))

    def fail(self, measRecord, error=None):
        measRecord.set(self.flag, True)
