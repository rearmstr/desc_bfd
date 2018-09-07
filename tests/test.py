from lsst.daf.persistence import Butler
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.geom as geom
import bfd
import galsim
import numpy as np

b=Butler('/datasets/hsc/repo/rerun/RC/w_2018_24/DM-14688/')
#exposure =b.get('deepCoadd_calexp',tract=9813,patch='4,4',filter='HSC-I')
#meas = b.get('deepCoadd_meas',tract=9813,patch='4,4',filter='HSC-I')

exposure = b.get('calexp',visit=1228, ccd=49)
meas = b.get('src',visit=1228, ccd=49)

ds=[]
for ii in range(1000):
    sns=[]
    sigmas=np.arange(0.2,2,0.1)
    for i in sigmas:
        weight = bfd.KSigmaWeight(4, i)

        measRecord = meas[ii]

        center = measRecord.getCentroid()

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
            continue
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
        jacobian2 = jacobian

        xref = ra -  (jacobian[0,0]*(measRecord.getX() - box.getMinX()) +
                      jacobian[0,1]*(measRecord.getY() - box.getMinY()))
        yref = dec - (jacobian[1,0]*(measRecord.getX() - box.getMinX()) +
                      jacobian[1,1]*(measRecord.getY() - box.getMinY()))


        sky_pos = exposure.getWcs().pixelToSky(measRecord.getCentroid())
        uvref = (sky_pos.getRa().asArcseconds(), sky_pos.getDec().asArcseconds())

        xref2 = ra -  (jacobian2[0,0]*(measRecord.getX() - box.getMinX()) +
                       jacobian2[0,1]*(measRecord.getY() - box.getMinY()))
        yref2 = dec - (jacobian2[1,0]*(measRecord.getX() - box.getMinX()) +
                       jacobian2[1,1]*(measRecord.getY() - box.getMinY()))
        local_x = measRecord.getX() - box.getMinX()
        local_y = measRecord.getY() - box.getMinY()

        bfd_wcs = bfd.WCS(jacobian, xyref=(local_x, local_y), uvref=uvref)
        bfd_wcs2 = bfd.WCS(jacobian2, xyref=(local_x, local_y), uvref=uvref)

        noise = measRecord.get("base_Variance_value")

        kdata = bfd.simpleImage(image, uvref, psf_image.array, wcs=bfd_wcs, pixel_noise=noise)
        kdata2 = bfd.simpleImage(image, uvref, psf_image.array, wcs=bfd_wcs2, pixel_noise=noise)
        moment_calc = bfd.MomentCalculator(kdata, weight, id=measRecord.getId())
        moment_calc2 = bfd.MomentCalculator(kdata2, weight, id=measRecord.getId())

        xyshift, error, msg = moment_calc.recenter()
        xyshift2, error2, msg2 = moment_calc2.recenter()
        moment = moment_calc.get_moment(0,0).all()
        moment2 = moment_calc2.get_moment(0,0).all()
        cov = moment_calc.get_covariance()[0]
        flux_cov = moment_calc.get_covariance()[0][0,0]
        sn_m = moment[0]/np.sqrt(flux_cov)
        sn_e = 0.5*(moment[2]/np.sqrt(cov[2,2])+(moment[3]/np.sqrt(cov[3,3])))
        sn_s =  moment[1]/np.sqrt(cov[1,1])
        sns.append(sn_s)

        
    if len(sns)==0:
        continue
    if np.all(np.isnan(sns)):
        continue
    best = np.nanargmax(sns)
    print (ii,sns[best],sigmas[best])
    ds.append(sigmas[np.argmax(sns)])
