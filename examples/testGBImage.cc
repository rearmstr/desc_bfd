#include "Image.h"
#include "Moments.h"
#include "PixelGalaxy.h"
#include "Interpolant.h"
#include "SBPixel.h"
#include "SbPsf.h"
#include "Image.h"
#include "Great3.h"
#include "KWeight.h"
#include "FitsImage.h"
using namespace std;
using namespace img;
using namespace old_bfd;
using namespace sbp;

int main(int argc,
         char *argv[])
{

    double flux = 100.;
    double sigma = 3.5;
    double e = 0.;
    double sigmaW = 1.;
    double beta = 0.;
    double x0 = 0.;
    double y0 = 0.;

    flux = 1.;
    sigma = 3.;
    e = 0.;//0.7;
    sigmaW = 0.5;
    beta = 0;//40. * 3.1415/180.;
    x0 = 0.;
    y0 = 0.;
    double noise = 1.;
    GaussianWeight kw(sigmaW);
    GaussianGalaxy<> gal(kw, flux, sigma*pow(1-e*e,-0.25), e, beta, x0, y0, noise);

    // Use delta-function PSF now.
    DeltaPsf psf;

    SBGaussian src(flux, sigma);
    Shear s;
    s.setEBeta(e,beta);
    Ellipse distortion(s, 0., Position<double>(x0,y0));
    SBProfile* lensed = src.distort(distortion);

    //**    double dx = 0.2;
    double dx = 1.;
  img::Image<> data = lensed->draw(dx);
//    img::Image<> data = src.draw(dx);
    data *= (dx*dx);    // Turn into SB in flux per pixel.
    cerr << "Image dimensions: " << data.getBounds() << endl;
    cerr << "Origin value: " << data(0,0) << endl;
    cout<<data.getBounds()<<endl;
    DMatrix22 a(0.);
    a(0,0) = a(1,1) = dx;
    DVector2 origin(0.);
    Affine map(a,origin);

    PixelGalaxy<15> pg(kw, data, map, psf, noise * (dx*dx));
    Moments<15> numeric = pg.getMoments();
    cout<<numeric<<endl;

    data.shift(1,1);
    img::FitsImage<>::writeToFITS("test.fits",data);
    

//     int max=15;
//     int maxPsf=25;
// //    flux=10;
//     Bounds<int> bounds(-max,max,-max,max);
//     Bounds<int> boundsPsf(0,maxPsf,0,maxPsf);
//     img::Image<> im(bounds);
//     img::Image<> imPsf(bounds);
//     double xc=0.5;
//     double yc=0.5;
//     double xcPsf=xc;
//     double ycPsf=yc;
//     //  sigma = 1.5;
//     double sigmaPsf = 1;
//     //sigmaW = 1.;
//     int ii=0;
//     for (int iy = im.yMin(); iy<= im.yMax(); iy++) {
//         for (int ix = im.xMin(); ix<= im.xMax(); ix++) {
//             double xx=ix;
//             double yy=iy;
//             im(ix,iy) = 1./(2*PI*sigma*sigma)*std::exp(-0.5*( (xx-xc)*(xx-xc)+(yy-yc)*(yy-yc))/(sigma*sigma));
//             ii+=1;
//         }
//     }
//     std::cout<<"Total "<<ii<<std::endl;
//     for (int iy = imPsf.yMin(); iy<= imPsf.yMax(); iy++) {
//         for (int ix = imPsf.xMin(); ix<= imPsf.xMax(); ix++) {
//             double xx=ix;
//             double yy=iy;
//             imPsf(ix,iy) = 1./(2*PI*sigmaPsf*sigmaPsf)*std::exp(-0.5*( (xx-xcPsf)*(xx-xcPsf)+(yy-ycPsf)*(yy-ycPsf))/(sigmaPsf*sigmaPsf));
//         }
//     }


//     fft::SincInterpolant interp1d;
//     fft::InterpolantXY interp2d(interp1d);

//     // Make an SBProfile from the stamp
//     sbp::SBPixel pixPsf(imPsf, interp2d);
//     pixPsf.setFlux(1.);

//     SbPsf sbPsf(pixPsf);
//     DVector2 origin2(0.);

//     origin2[0] = xc;
//     origin2[1] = yc;
//     Affine map2(a,origin2);

//     //im.shift(-max,-max);
//     //cout<<im(0,0)<<endl;
//     PixelGalaxy<> pg2(kw, im, map2, psf,noise);
//     cout<<pg2.getMoments()<<endl;
//     imPsf.shift(1,1);
//     img::FitsImage<>::writeToFITS("test.fits",imPsf);
//     cout<<pg2.getCov()<<endl;

//     vector<double> kvals(25);
//     for(int i=0;i<25;++i) kvals[i]=i*0.1777;
//     vector<double> p(25);
//     p[0]=0.000000;
//     p[1]=1.719067;
//     p[2]=1.433970;
//     p[3]=1.359070;
//     p[4]=1.271386;
//     p[5]=1.233483;
//     p[6]=1.233999;
//     p[7]=1.229341;
//     p[8]=1.221745;
//     p[9]=1.233884;
//     p[10]=1.205419;
//     p[11]=1.167032;
//     p[12]=1.141703;
//     p[13]=1.068160;
//     p[14]=0.997798;
//     p[15]=0.924048;
//     p[16]=0.850920;
//     p[17]=0.781266;
//     p[18]=0.763041;
//     p[19]=0.717132;
//     p[20]=0.655426;
//     p[21]=0.573232;
//     p[22]=0.500142;
//     p[23]=0.442679;
//     p[24]=0.409335;

    //Table<double,double> table(kvals, p);
    //cout<<pg2.getCov(&table)<<endl;
}

// }
