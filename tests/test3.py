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

exposure =b.get('calexp',visit=1228, ccd=49)
meas = b.get('src',visit=1228, ccd=49)


moms=[]
for i in range(6,7):
    weight = bfd.KSigmaWeight(4, 13)
    print(i)
    measRecord = meas[i]
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
        moms.append(np.zeros(6))
        continue

    image = exposure.image[box].array

    # PSF image must also be the same size
    psf_image_base = galsim.ImageF(exposure.getPsf().computeKernelImage(center).array)
    psf_image_interp = galsim.InterpolatedImage(psf_image_base, scale=1)
    psf_image = galsim.ImageF(box_size, box_size)
    psf_image_interp.drawImage(psf_image, method='no_pixel')

    ra = box.getMinX()#measRecord.get('coord_ra').asArcseconds()
    dec = box.getMinY()#measRecord.get('coord_dec').asArcseconds()        

    local_lin_wcs = exposure.getWcs().linearizePixelToSky(center, geom.arcseconds)
    #jacobian = local_lin_wcs.getLinear().getMatrix()
    jacobian = np.array([[1,0],[0,1]])
    xref = ra -  (jacobian[0,0]*(measRecord.getX() - box.getMinX()) +
                  jacobian[0,1]*(measRecord.getY() - box.getMinY()))
    yref = dec - (jacobian[1,0]*(measRecord.getX() - box.getMinX()) +
                  jacobian[1,1]*(measRecord.getY() - box.getMinY()))


    sky_pos = exposure.getWcs().pixelToSky(measRecord.getCentroid())
    uvref = (xref,yref)#(sky_pos.getRa().asArcseconds(), sky_pos.getDec().asArcseconds())

    bfd_wcs = bfd.WCS(jacobian, xyref=(0,0), uvref=(0,0))

    noise = 1#measRecord.get("base_Variance")

    kdata = bfd.simpleImage(image, (15.5,15.5), psf_image.array, wcs=bfd_wcs, pixel_noise=noise)
    
    moment_calc = bfd.MomentCalculator(kdata, weight, id=measRecord.getId())

    moment = moment_calc.get_moment(0,0).all()

    xyshift, error, msg = moment_calc.recenter()
    moment = moment_calc.get_moment(0,0).all()
    moms.append(moment)


 
