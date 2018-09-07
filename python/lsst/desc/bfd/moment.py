
import lsst.pex.config as pexConfig
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.geom as geom
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin
from lsst.meas.base.pluginRegistry import register
import numpy as np

import bfd
import galsim

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
        self.cov_even = schema.addField('bfd_cov_even', type="ArrayF", size=NMoment*(NMoment+1)//2, doc="even moment covariance matrix")
        self.cov_odd = schema.addField('bfd_cov_odd', type="ArrayF", size=3, doc="odd moment covariance matrix")
        self.flag = schema.addField('bfd_flag', type="Flag", doc="Set to 1 for any fatal failure")
        self.centroid_flag = schema.addField('bfd_flag_centroid', type="Flag", doc="Set to 1 for any fatal failure of centroid")

        self.weight = bfd.KSigmaWeight(self.config.weight_n, self.config.weight_sigma)


    def measure(self, measRecord, exposure):
        center = measRecord.getCentroid()

        # currently box must be square for python based measurements
        orig_box = measRecord.getFootprint().getBBox()
        w, h = orig_box.getWidth(), orig_box.getHeight()
        box_size = max(w, h)
        if box_size % 2 == 1:
            box_size += 1
        box = afwGeom.Box2I(afwGeom.Point2I(center.getX()-box_size/2, center.getY()-box_size/2), afwGeom.Extent2I(box_size, box_size))

        if not exposure.getBBox().contains(box):
            box_size = min(w, h)
            if box_size % 2 == 1:
                box_size += 1
            box = afwGeom.Box2I(afwGeom.Point2I(center.getX()-box_size/2, center.getY()-box_size/2), afwGeom.Extent2I(box_size, box_size))

        if not exposure.getBBox().contains(box):
            measRecord.set(self.flag, 1)
            return

        image = exposure.image[box].array

        # PSF image must also be the same size
        psf_image_base = galsim.ImageF(exposure.getPsf().computeKernelImage(center).array)
        psf_image_interp = galsim.InterpolatedImage(psf_image_base, scale=1)
        psf_image = galsim.ImageF(box_size, box_size)
        psf_image_interp.drawImage(psf_image, method='no_pixel')

        ra = measRecord.get('coord_ra').asArcseconds()
        dec = measRecord.get('coord_dec').asArcseconds() 
        local_lin_wcs = exposure.getWcs().linearizePixelToSky(center, geom.arcseconds)

        jacobian = local_lin_wcs.getLinear().getMatrix()
        sky_pos = exposure.getWcs().pixelToSky(center)
        uvref = (sky_pos.getRa().asArcseconds(), sky_pos.getDec().asArcseconds())

        xy_pos = (center.getX() - box.getMinX(), center.getY() - box.getMinY())
        bfd_wcs = bfd.WCS(jacobian, xyref=xy_pos, uvref=uvref)
        
        noise = measRecord.get("base_Variance_value")

        kdata = bfd.simpleImage(image, uvref, psf_image.array, wcs=bfd_wcs, pixel_noise=noise)
        moment_calc = bfd.MomentCalculator(kdata, self.weight, id=measRecord.getId())
        
        xyshift, error, msg = moment_calc.recenter()
        if error:
            measRecord.set(self.flag, 1)
            measRecord.set(self.centroid_flag, 1)
            return
        else:
            #cov_even, cov_odd = moment_calc.get_covariance()
            covgal = moment_calc.get_covariance()
            moment = moment_calc.get_moment(0,0)
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
        #self.log.info('BFD failed', error)
        #print(error)
        measRecord.set(self.flag, True)
        assert(False)
