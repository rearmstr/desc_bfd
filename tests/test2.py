import galsim
import bfd
import numpy as np



def simpleImage(image, origin, psf, pixel_scale=1.0, pixel_noise=None, wcs=None, noise_ps=None):
    '''Create PSF-corrected k-space image and variance array
    image  postage stamp to use, 2d numpy array, in units of FLUX PER PIXEL.
           2nd axis taken as x.
    origin Location of the origin of the galaxy (2 element array).  Should be near
           the image center. If a wcs is given, these are in the world coordinates
           of the WCS.  If wcs is None, these are 0-index pixel coordinates.
    psf    matching-size postage stamp for the PSF.  We will fix normalization.
           Assumed to have origin at [N/2,N/2] pixel.
    pixel_scale is number of sky units per pixel, only used if wcs is None
    pixel_noise RMS noise of the image pixel values.  
    wcs    A WCS instance giving the map from 0-indexed pixel coordinates to a world
           (sky) system.

    Returns KData instance.  All quantities therein are in sky units and the flux units of
           the input image.
    '''

    if psf.shape != image.shape or psf.ndim != 2 or image.shape[0]!=image.shape[1]:
        raise Exception('PSF and image shape must be matching square arrays for simpleImage')
    # ??? Do zero padding to some nominal size(s)
    N = image.shape[0]
    if not N%2==0:
        raise Exception('Image and PSF need to be even-dimensioned.')
    kval =  np.fft.rfft2(image)
    # Flip signs to shift coord origin from [0,0] to N/2,N/2
    kval[1::2,::2] *= -1.
    kval[::2,1::2] *=-1.
    
    # Same for PSF
    kpsf  =  np.fft.rfft2(psf)
    kpsf[1::2,::2] *= -1.
    kpsf[::2,1::2] *=-1.
    # Normalize PSF to unit flux (note DC is at [0,0])
    kpsf /= kpsf[0,0]

    # Correct for PSF
    kval /= kpsf
    
    # Double the weight on samples whose conjugates are missing
    conjugate = np.zeros_like(kpsf, dtype=bool)
    conjugate[:,1:N//2] = True

    ky = np.ones_like(kpsf, dtype=float) \
      * np.fft.fftfreq(N)[:,np.newaxis]*2.0*np.pi
    kx = np.ones_like(kpsf, dtype=float) \
      * np.fft.rfftfreq(N)*2.0*np.pi

    # K area per sample
    d2k = (kx[0,1] - kx[0,0])**2


        # Make variance array if we have noise
    kvar = None
    if pixel_noise is not None and noise_ps is None:
        kvar = np.ones_like(kpsf,dtype=float) * (pixel_noise*N)**2
        kvar /= (kpsf.real*kpsf.real + kpsf.imag*kpsf.imag)
    elif pixel_noise is None and noise_ps is not None:
        noise = np.sqrt(noise_ps(np.sqrt(kx**2 + ky**2)))
        kvar = np.ones_like(kpsf,dtype=float) * (noise*N)**2
        kvar /= (kpsf.real*kpsf.real + kpsf.imag*kpsf.imag)


    # Adjust kx, ky, d2k for coordinate mapping
    # and calculate (sky coordinate) displacement of
    # galaxy origin from FFT phase center
    if wcs is None:
        dxy = np.array(origin) - N/2  # origin was given in pixels
        kx /= pixel_scale
        ky /= pixel_scale
        d2k /= pixel_scale**2
    else:
        # Give needed origin shift in sky units
        dxy = np.array(origin) - wcs.getuv( np.array([-N/2, -N/2],dtype=float))
        # Put k values into an array of shape (2, N^2)
        kxy = np.vstack((kx.flatten(),ky.flatten()))
        # Transform k's and put back into shape of the k arrays
        kxy = wcs.getuv_k(kxy)
        kx = kxy[0].reshape(kx.shape)
        ky = kxy[1].reshape(ky.shape)
        
        # Rescale pixel area in k space, need absolute value of determinant 
        # to keep d2k positive
        d2k /= np.abs(wcs.getdet())

    # Apply phase shift to move center
    phase = kx * dxy[0] + ky * dxy[1]
    
    kval *=np.exp(1j*phase)

    return kval,kx,ky,d2k,conjugate,kvar

def rshift(array):
    N = array.shape[0]
    shift = np.zeros_like(array)
    shift[N//2:]=array[:N//2]
    shift[:N//2]=array[N//2:]

    return shift



gauss = galsim.Gaussian(sigma=5,flux=300)
gauss = gauss.shear(g1=0.05,g2=-0.1)
psf = galsim.Gaussian(sigma=0.5,flux=1)
gal = galsim.Convolve(psf, gauss)
N=100
scale=1
d11=-scale
d22=-scale
d12=0.1
d21=-0.5
sky_center = galsim.CelestialCoord(ra=0*galsim.degrees, dec=0*galsim.degrees)
affine1 = galsim.AffineTransform(d11,d12,d21,d22, origin=galsim.PositionD((N/2),(N/2)))

wcs1 = galsim.TanWCS(affine1, sky_center)


image1 = galsim.ImageF(N,N)
psf_image1 = galsim.ImageF(N,N)


psf.drawImage(image=psf_image1, wcs=wcs1)

gal.drawImage(image=image1, wcs=wcs1)

weight = bfd.KSigmaWeight(4, 2*scale)
uvref=(0,0)
jacobian1 = np.array([[d11,d12],[d21,d22]])

bfd_wcs = bfd.WCS(jacobian1, xyref=(N/2, N/2), uvref=uvref)

noise=1
kval,kx,ky,d2k,conjugate,kvar = simpleImage(image1.array, uvref , psf_image1.array, wcs=bfd_wcs)
kshift=rshift(kval)
weight.set_k(kx,ky)
wshift=rshift(weight.w)

kdata2 = bfd.simpleImage(image1.array, uvref, psf_image1.array, wcs=bfd_wcs, pixel_noise=noise)
moment_calc2 = bfd.MomentCalculator(kdata2, weight)

moment2 = moment_calc2.get_moment(0,0).all()
print(moment2)

