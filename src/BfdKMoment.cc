// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2015 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */


#include "lsst/desc/bfd/BfdKMoment.h"
#include "lsst/desc/bfd/ImageConvert.h"

//#include "lsst/log/Log.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/daf/base.h"

#include "Moments.h"
#include "PixelGalaxy.h"
#include "Interpolant.h"
#include "SBPixel.h"
#include "Prior.h"
#include "KdPrior.h"
#include "SbPsf.h"
#include "Image.h"
#include "Great3.h"
#include "KWeight.h"
#include "FitsImage.h"
#include "Random.h"
//#include "MultiGalaxy.h"

#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

// Rename namespace to avoid collision
namespace BFD = old_bfd;


namespace lsst { namespace desc { namespace bfd {


 ////////////////////////////////////////////////////////////////////////////////////////////
 //-------------------- Utility Code --------------------------------------------------------
 ////////////////////////////////////////////////////////////////////////////////////////////

 namespace {
 ran::UniformDeviate myrand;
 PTR(BFD::PixelGalaxy<UseMoments>)
 getPixelGalaxy(
     BFD::KWeight &kw,
     PTR(afw::image::Image<Pixel>) const & image,
     CONST_PTR(afw::detection::Psf) const & psf,
     afw::geom::LinearTransform const & transform,
     afw::geom::Point2D const & center,
     double noise,
     double &x0,
     double &y0,
     bool reCentroid,
     bool reCentroidPsf=true,
     bool ignorePsf=false,
     int interpOrder=5,
     //Table<> noisePs=Table<>(),
     vector<double> kval=vector<double>(),
     vector<double> val=vector<double>(),
     PTR(afw::image::Image<Pixel>) cov=PTR(afw::image::Image<Pixel>)()
     )
 {
     LOG_LOGGER trace3Logger = LOG_GET("TRACE3.desc.bfd");

     ImageConverter<Pixel> imageConvert(image);
     LOGL_DEBUG(trace3Logger, "image size %d %d",image->getBBox().getWidth(),image->getBBox().getHeight() );
     PTR(img::Image<Pixel>) im = imageConvert.getBFDImage();
     DVector2 origin(0.);
     if (reCentroid)
       {
         double hlr;
         double sigma = image->getWidth() / 12.;
         BFD::centroid(*im, sigma, x0, y0, hlr);
         LOGL_DEBUG(trace3Logger,"Calculated centroid %g %g, half light radius %g",x0,y0,hlr);
     } else {
       LOGL_DEBUG(trace3Logger,"Using centroid %g %g",x0,y0);
     }

     origin[0] = x0;
     origin[1] = y0;

     // Create the affine transform
     DMatrix22 a;
     Eigen::Matrix2d in = transform.getMatrix();
     std::copy(in.data(), in.data() + 4, a.ptr());
     BFD::Affine map(a,origin);
     double detA = std::abs(map.jacobian());

     double scaling = std::sqrt(std::abs(detA));

     fft::Lanczos interp1d(interpOrder,true);
     fft::InterpolantXY interp2d(interp1d);
     // Useful for debugging images
     ////if (debugLog.getThreshold()<-5) {
     ////        stringstream name;
     ////        int id = (10000*myrand());
     ////        name<<"test_image"<<id<<".fits"<<ends;
     ////        image->writeFits(name.str());
     ////        stringstream psf_name;
     ////        psf_name<<"test_psf_image"<<id<<".fits"<<ends;
     ////        psf->computeKernelImage(afw::geom::Point2D(x0,y0))->writeFits(psf_name.str());
     ////}

     PTR(BFD::Psf) sbPsf;
      if (ignorePsf) {
          sbPsf.reset(new BFD::DeltaPsf());
      } else {
          // The image defined in psf->computeImage depends on what is defined in afw::math::Kernel
          PTR(afw::image::Image<Pixel>) psfImage(
                      new afw::image::Image<Pixel>(*psf->computeKernelImage(afw::geom::Point2D(x0,y0)), true));

          ImageConverter<Pixel> imageConvertPsf(psfImage);
          PTR(img::Image<Pixel>) psfIm = imageConvertPsf.getBFDImage();
          // Make an SBProfile from the stamp
          sbp::SBPixel pixPsf(*psfIm, interp2d);
          pixPsf.setFlux(1.);
          sbPsf.reset(new BFD::SbPsf(pixPsf));

          if (reCentroidPsf) {
 	   ///LOGL_DEBUG(trace3Logger,"Recentroiding PSF");
              double psfHlr, psfX0(0.), psfY0(0.);
              double psfSigma = psfImage->getWidth() / 12.;
              BFD::centroid(*psfIm, psfSigma, psfX0, psfY0, psfHlr);
              LOGL_DEBUG(trace3Logger,"Calculated psf centroid %g %g, half light radius %g",psfX0,psfY0,psfHlr);

              int xOriginFT = (psfIm->xMin() + psfIm->xMax() + 1)/2;
              int yOriginFT = (psfIm->yMin() + psfIm->yMax() + 1)/2;

              // Change FT phases to move centroid back to the FFT phase center
              LOGL_DEBUG(trace3Logger,"Shift psf %g %g",xOriginFT-psfX0, yOriginFT-psfY0);
              sbp::SBProfile* shifted = pixPsf.shift(xOriginFT-psfX0, yOriginFT-psfY0);
              BFD::SbPsf out(*shifted);
              sbPsf.reset(new BFD::SbPsf(out));
              delete shifted;
          }
      }

      LOGL_DEBUG(trace3Logger,"PSF kMax %g",sbPsf->kMax());
      LOGL_DEBUG(trace3Logger,"Weight kMax %g",kw.kMax());

     // For now avoid doing wcs parts
     //(*im) *= detA;
     //noise *= detA;
     if (kval.size() > 0) {
       LOGL_DEBUG(trace3Logger,"getting PixelGalaxy using noise power spectrum");
         PTR(BFD::PixelGalaxy<UseMoments>) pg(new BFD::PixelGalaxy<UseMoments>(kw, *im, map, *sbPsf, noise));
         pg->setNoisePs(kval, val);
         return pg;
     } else if (cov) {
       LOGL_DEBUG(trace3Logger,"getting PixelGalaxy using noise power image");

       assert(false);
         ImageConverter<Pixel> imageConvertCov(cov);
         PTR(img::Image<Pixel>) covIm = imageConvertCov.getBFDImage();

         // Make an SBProfile from the stamp
         sbp::SBPixel covPix(*covIm, interp2d);
         sbp::SBProfile* shifted = covPix.shift(-9, -9);

         //PTR(BFD::PixelGalaxy<UseMoments>) pg(new BFD::PixelGalaxy<UseMoments>(kw, *im, map, *sbPsf, &covPix, noise));
         //return pg;
     }

     PTR(BFD::PixelGalaxy<UseMoments>) pg(new BFD::PixelGalaxy<UseMoments>(kw, *im, map, *sbPsf, noise));
     return pg;
 }
 }




 ////////////////////////////////////////////////////////////////////////////////////////////
 //-------------------- Bfd KMoment Code --------------------------------------------------------
 ////////////////////////////////////////////////////////////////////////////////////////////
 BfdKMomentResultKey BfdKMomentResultKey::addFields(
     afw::table::Schema & schema,
     std::string const & name
     )
 {
     BfdKMomentResultKey r;
     r._momentsKey = schema.addField<afw::table::Array<float> >(name+".moments", "k-space moments", NMoment);
     r._momentsCovKey = schema.addField<afw::table::Array<float> >(name+".momentsCov",
                                                                   "covariance of k-space moments",
                                                                   NMoment*NMoment);
     r._momentsPsfKey = schema.addField<afw::table::Array<float> >(name+".moments.psf",
                                                                   "k-space moments of Psf", NMoment);
     r._centerKey = afw::table::PointKey<double>::addFields(schema, name+".center",
 							   "center where moments were evaluated", "pixel");
     r._shiftKey = afw::table::PointKey<double>::addFields(schema, name+".shift",
 							  "amount center was shifted", "pixel");
     r._flagsKey[0] = schema.addField<afw::table::Flag>(name+".flags", "failure flag");
     r._flagsKey[1] = schema.addField<afw::table::Flag>(name+".flags.shift.failed", "shift failed");
     r._flagsKey[2] = schema.addField<afw::table::Flag>(name+".flags.shift.large", "shift too large");
     r._flagsKey[3] = schema.addField<afw::table::Flag>(
         name + ".flags.shift.centroid.large",
         "centroid of shifted moments too large");
     r._flagsKey[4] = schema.addField<afw::table::Flag>(name+".flags.too-big", "pixel region too big");
     r._flagsKey[5] = schema.addField<afw::table::Flag>(name+".flags.variance",
                                                        "problem with variance calculation");
     r._flagsKey[6] = schema.addField<afw::table::Flag>(name+".flags.parent",
                                                        "parent with children");
     r._flagsKey[7] = schema.addField<afw::table::Flag>(name+".flags.flux_negative",
                                                        "bad moment flux");
     r._flagsKey[8] = schema.addField<afw::table::Flag>(name+".flags.size_negative",
                                                        "bad moment size");
     r._flagsKey[9] = schema.addField<afw::table::Flag>(name+".flags.saturated.center",
                                                        "saturated center");
     r._flagsKey[10] = schema.addField<afw::table::Flag>(name+".flags.footprint-empty",
                                                        "empty footprint");
     return r;

 }

 BfdKMomentResultKey::BfdKMomentResultKey(afw::table::Schema schema, std::string name) {
     _momentsKey = schema.find<afw::table::Array<float> >(name+".moments").key;
     _momentsCovKey = schema.find<afw::table::Array<float> >(name+".momentsCov").key;
     _momentsPsfKey = schema.find<afw::table::Array<float> >(name+".moments.psf").key;
     _centerKey = lsst::afw::table::PointKey<double>(schema[name+".center"]);
     _shiftKey = lsst::afw::table::PointKey<double>(schema[name+".shift"]);
     _flagsKey[0] = schema.find<afw::table::Flag>(name+".flags").key;
     _flagsKey[1] = schema.find<afw::table::Flag>(name+".flags.shift.failed").key;
     _flagsKey[2] = schema.find<afw::table::Flag>(name+".flags.shift.large").key;
     _flagsKey[3] = schema.find<afw::table::Flag>(name+".flags.shift.centroid.large").key;
     _flagsKey[4] = schema.find<afw::table::Flag>(name+".flags.too-big").key;
     _flagsKey[5] = schema.find<afw::table::Flag>(name+".flags.variance").key;
     _flagsKey[6] = schema.find<afw::table::Flag>(name+".flags.parent").key;
     _flagsKey[7] = schema.find<afw::table::Flag>(name+".flags.flux_negative").key;
      _flagsKey[8] = schema.find<afw::table::Flag>(name+".flags.size_negative").key;
      _flagsKey[9] = schema.find<afw::table::Flag>(name+".flags.saturated.center").key;
     _flagsKey[10] = schema.find<afw::table::Flag>(name+".flags.footprint-empty").key;
 }

 void BfdKMomentResultKey::set(afw::table::BaseRecord & record, BfdKMomentResult const & value) const {
     record.set(_momentsKey, value.moments);
     record.set(_momentsCovKey, ndarray::flatten<1>(value.momentsCov));
     record.set(_momentsPsfKey, value.momentsPsf);
     _centerKey.set(record, value.center);
     _shiftKey.set(record, value.shift);
     // for (int b = 0; b < BfdKMomentResult::N_FLAGS; ++b) {
     //     record.set(_flagsKey[b], value.flags[b]);
     // }
 }

 BfdKMomentResult BfdKMomentResultKey::get(afw::table::BaseRecord const & record) const {
     BfdKMomentResult result;
     ndarray::flatten<1>(result.moments) = record.get(_momentsKey);
     ndarray::flatten<1>(result.momentsCov) = record.get(_momentsCovKey);
     ndarray::flatten<1>(result.momentsPsf) = record.get(_momentsPsfKey);
     result.center = _centerKey.get(record);
     result.shift = _shiftKey.get(record);
     result.shift = record.get(_shiftKey);
     // Do not add flags for now because we need something better
     // for (int b = 0; b < BfdKMomentResult::N_FLAGS; ++b) {
     //   result.flags[b] = record.get(_flagsKey[b]);
     // }

     return result;
 }

 BfdKMoment::BfdKMoment(
     Control const & ctrl,
     std::string const & name,
     afw::table::Schema & schema
 ) :
     _ctrl(ctrl),
     _resultKey(ResultKey::addFields(schema, "bfd"))
 {

     if ( !ctrl.calculateVariance && !ctrl.useTableVariance && !ctrl.useRecVariance ) {
         throw LSST_EXCEPT(
             pex::exceptions::InvalidParameterError,
             "Must set either calculateVariance, useTableVariance or useRecVariance"
         );
     }

     if ( ctrl.calculateVariance + ctrl.useTableVariance + ctrl.useRecVariance > 1 ) {
         throw LSST_EXCEPT(
             pex::exceptions::InvalidParameterError,
             "Can only set one from calculateVariance,useTableVariance or useRecVariance"
         );
     }
 }

 BfdKMoment::Result BfdKMoment::measure_image(
     BfdKMomentControl const & ctrl,
     PTR(afw::image::Image<Pixel>)const & image,
     CONST_PTR(afw::detection::Psf) const & psf,
     afw::geom::LinearTransform const & transform,
     afw::geom::Point2D const & center,
     double noise
 ) {

     BfdKMoment::Result result;
     LOG_LOGGER trace3Logger = LOG_GET("TRACE3.desc.bfd");
     LOG_LOGGER trace4Logger = LOG_GET("TRACE4.desc.bfd");

     double scaling = std::sqrt(std::abs(transform.computeDeterminant()));
     PTR(BFD::KWeight) kw;
     if (ctrl.wIndex < 0) {
         kw.reset(new BFD::GaussianWeight(ctrl.sigma*scaling));
     } else {
         kw.reset(new BFD::KSigmaWeight(ctrl.sigma*scaling, ctrl.wIndex));
     }

     try {
         result.flags[0] = false;
         double x0=center.getX();
         double y0=center.getY();


         PTR(BFD::PixelGalaxy<UseMoments>) pg = getPixelGalaxy(*kw, image, psf, transform, center, noise, x0, y0,
                                                               ctrl.reCentroid, ctrl.reCentroidPsf, ctrl.ignorePsf, ctrl.interpOrder);
         BFD::Moments<UseMoments> moments;
         BFD::MomentCovariance<UseMoments> momentsCov;

         DVector2 origin(0.);
         origin[0] = x0;
         origin[1] = y0;
         moments = pg->getMoments();
         momentsCov = pg->getCov();


         LOGL_DEBUG(trace3Logger,"Using center %g %g", origin[0], origin[1]);
         LOGL_DEBUG(trace4Logger,"got moments: %g", moments[0]);
         LOGL_DEBUG(trace4Logger,"             %g", moments[1]);
         LOGL_DEBUG(trace4Logger,"             %g", moments[2]);
         LOGL_DEBUG(trace4Logger,"             %g", moments[3]);
         LOGL_DEBUG(trace4Logger,"             %g", moments[4]);
         LOGL_DEBUG(trace4Logger,"             %g", moments[5]);
         LOGL_DEBUG(trace4Logger,"flux S/N: %g",moments[0]/std::sqrt(momentsCov(0,0)));
         LOGL_DEBUG(trace4Logger,"cov(0,0): %g",momentsCov(0,0));

         BFD::Moments<UseMoments> momentsPsf;
         if (!ctrl.ignorePsf) {
             LOGL_DEBUG(trace3Logger,"Get moments of the psf");
             // Get the moments of the psf.  Use computeImage so that the psf will be centered on
             // the point x0,y0
             PTR(afw::image::Image<Pixel>) psfImage(
                 new afw::image::Image<Pixel>(*psf->computeImage(afw::geom::Point2D(x0,y0)), true));
             bool reCentroid = false;
             bool reCentroidPsf = false;
             bool ignorePsf = true;
             PTR(BFD::PixelGalaxy<UseMoments>) pgPsf = getPixelGalaxy(*kw, psfImage, psf, transform,
                                                                      center, 1, x0, y0,
                                                                      reCentroid, reCentroidPsf, ignorePsf, ctrl.interpOrder);

             momentsPsf=pgPsf->getMoments();
             LOGL_DEBUG(trace4Logger,"Using center for PSF %g %g", x0, y0);
             LOGL_DEBUG(trace4Logger,"got moments for the PSF: %g", momentsPsf[0]);
             LOGL_DEBUG(trace4Logger,"                         %g", momentsPsf[1]);
             LOGL_DEBUG(trace4Logger,"                         %g", momentsPsf[2]);
             LOGL_DEBUG(trace4Logger,"                         %g", momentsPsf[3]);
             LOGL_DEBUG(trace4Logger,"                         %g", momentsPsf[4]);
             LOGL_DEBUG(trace4Logger,"                         %g", momentsPsf[5]);
         }


         if (ctrl.shift) {
             double dx, dy;
             PTR(BFD::TemplateGalaxy<UseMoments>) shifted;
             try {
                 shifted.reset(BFD::newtonShift(*pg, dx, dy, ctrl.shiftIterations));
             } catch(std::runtime_error const &error) {
                 LOGL_DEBUG(trace3Logger,"Shift failed");
                 result.flags[BfdKMomentResult::SHIFT_FAILED] = true;
                 throw;
             }

             if ( (dx*dx + dy*dy) > ctrl.maxShift*ctrl.maxShift) {
                 LOGL_DEBUG(trace3Logger,"Cenroid shift to null moments too large %g", std::sqrt(dx*dx + dy*dy) );
                 result.flags[BfdKMomentResult::SHIFT_LARGE] = true;
                 //result.flags[0] = true;
                 //return result;
             }
             LOGL_DEBUG(trace3Logger,"Shift to null centroid moments %g %g",dx,dy);
             origin[0] += dx;
             origin[1] += dy;
             moments = shifted->getMoments();
             momentsCov = shifted->getCov();

             LOGL_DEBUG(trace4Logger,"new moments: %g",moments[0]);
             LOGL_DEBUG(trace4Logger,"             %g",moments[1]);
             LOGL_DEBUG(trace4Logger,"             %g",moments[2]);
             LOGL_DEBUG(trace4Logger,"             %g",moments[3]);
             LOGL_DEBUG(trace4Logger,"             %g",moments[4]);
             LOGL_DEBUG(trace4Logger,"             %g",moments[5]);
             LOGL_DEBUG(trace4Logger,"flux S/N: %g",moments[0]/std::sqrt(momentsCov(0,0)));
             LOGL_DEBUG(trace4Logger,"cov(0,0): %g",momentsCov(0,0));

             if ( std::abs(moments[1]) > ctrl.maxCentroid ||
                  std::abs(moments[2]) > ctrl.maxCentroid ) {
                 LOGL_DEBUG(trace3Logger,"Shifted centroid moments too large %g %g, flagging",moments[1],moments[2]);
                 result.flags[BfdKMomentResult::SHIFT_CENTROID_LARGE] = true;
                 //result.flags[0] = true;
             }

             result.shift.setX(dx);
             result.shift.setY(dy);

         }
         result.center.setX(origin[0]);
         result.center.setY(origin[1]);

         if ( moments[0] < 0 ) {
             LOGL_DEBUG(trace3Logger,"Flux < 0, flagging");
             result.flags[BfdKMomentResult::FLUX_NEGATIVE] = true;
             result.flags[0] = true;
         }

         if ( moments[3] < 0 ) {
             LOGL_DEBUG(trace3Logger,"Size < 0, flagging");
             result.flags[BfdKMomentResult::SIZE_NEGATIVE] = true;
             result.flags[0] = true;
         }

         std::copy(moments.ptr(), moments.ptr()+NMoment, result.moments.getData());
         std::copy(momentsCov.ptr(), momentsCov.ptr()+NMoment*NMoment, result.momentsCov.getData());
         if (!ctrl.ignorePsf) {
             std::copy(momentsPsf.ptr(), momentsPsf.ptr()+NMoment, result.momentsPsf.getData());
         }


     }  catch(std::runtime_error const &error) {
         LOGL_DEBUG(trace3Logger,"Runtime error: %s", error.what());
         result.flags[0] = true;
     } catch(...) {
         LOGL_DEBUG(trace3Logger,"Unknown exception, moment failed");
         result.flags[0] = true;
     }

     return result;
 }

 BfdKMoment::Result BfdKMoment::measure_exp(
     BfdKMomentControl const & ctrl,
     afw::geom::Box2I box,
     afw::image::Exposure<Pixel> const & exposure,
     afw::geom::Point2D const & center,
     double variance
 ) {

     PTR(afw::detection::Psf const) psf = exposure.getPsf();
     PTR(afw::image::Image<Pixel>) imageExp = exposure.getMaskedImage().getImage();
     PTR(afw::image::Image<Pixel>) image =
         std::make_shared<afw::image::Image<Pixel> >(*imageExp, box,
                                                       afw::image::PARENT);
     CONST_PTR(afw::geom::SkyWcs) wcs = exposure.getWcs(); 
     if (false) {
       //afw::geom::LinearTransform transform = wcs->linearizePixelToSky(center).getLinear();
       //lsst::geom::SpherePoint skyCoord = wcs->pixelToSky(center);
       //afw::geom::Point2D sky(skyCoord->toIcrs().getRa().asDegrees(), skyCoord->toIcrs().getDec().asDegrees());
     }

     afw::geom::LinearTransform transform=afw::geom::LinearTransform::makeScaling(1.);
     return measure_image(ctrl, image, psf, transform, center, variance);

 }


 void BfdKMoment::measure(
 	afw::table::SourceRecord & measRecord,
         afw::image::Exposure<Pixel> const & exposure) const
 {

   Control const& control = static_cast <Control const&> ( _ctrl);
   afw::geom::Point2D const & center = measRecord.getCentroid();
   this->_apply(measRecord, exposure, center);

 }


 void BfdKMoment::fail(
     afw::table::SourceRecord & record,
     meas::base::MeasurementError * error) const
 {
   BfdKMoment::Result result;
   result.flags[0] = true;
   _resultKey.set(record, result);
   return;
 }


 void BfdKMoment::_apply(
     afw::table::SourceRecord & source,
     afw::image::Exposure<Pixel> const & exposure,
     afw::geom::Point2D const & center
 ) const {

     Control const& control = static_cast <Control const&> ( _ctrl);
     LOG_LOGGER trace3Logger = LOG_GET("TRACE3.desc.bfd");
     LOG_LOGGER trace4Logger = LOG_GET("TRACE4.desc.bfd");

     double variance = 1;
     std::vector<double> kData, psData;

     // Calculate the variance in the image
     try {
         if (control.calculateVariance) {
             // Cache the noise for different images.  There is not a good way to pass this in
             // if we use the measurement framework.  Think about building a singleton to handle this
             if (_noiseCache.find(exposure.getId()) == _noiseCache.end()) {
                 LOGL_DEBUG(trace3Logger,"Current id not in cache");
                 afw::math::StatisticsControl sctrl;
                 sctrl.setNumIter(3);
                 sctrl.setNumSigmaClip(5.0);
                 sctrl.setAndMask(0x0);
                 sctrl.setNanSafe(true);
                 afw::math::Statistics stats =
                     afw::math::makeStatistics(*exposure.getMaskedImage().getImage(), afw::math::VARIANCECLIP);
                 double var = stats.getValue(afw::math::VARIANCECLIP);
                 LOGL_DEBUG(trace3Logger,"Calculated variance is %g for id=%d", var, exposure.getId());
                 _noiseCache[exposure.getId()] = var;

             }
             LOGL_DEBUG(trace3Logger,"Using cache to calculate variance");
             variance = _noiseCache[exposure.getId()];
         }
         else if (control.useRecVariance) {
             LOGL_DEBUG(trace3Logger,"Using variance from catalog");
             afw::table::Key<double> noiseKey = source.getSchema()["noise_variance"];
             variance =  source.get(noiseKey);
             LOGL_DEBUG(trace3Logger,"Using variance from catalog %g",variance);
         }
         else if (control.useTableVariance) {
             PTR(daf::base::PropertyList) meta = source.getTable()->getMetadata();
             variance = meta->get<double>("noise_variance");
             LOGL_DEBUG(trace3Logger,"Using variance from meta data %g",variance);
         }

     } catch(const pex::exceptions::Exception &ex) {
    
         throw;
     }
     LOGL_DEBUG(trace3Logger,"Variance is %g %d",variance,control.useNoisePs);

     // setup images and psf, I am not currently using the wcs calculation
     CONST_PTR(afw::detection::Psf) psf = exposure.getPsf();
     PTR(afw::image::Image<Pixel>) imageExp = exposure.getMaskedImage().getImage();
     PTR(afw::image::Image<Pixel>) image =
         std::make_shared<afw::image::Image<Pixel> >(*imageExp, source.getFootprint()->getBBox(),
                                                       afw::image::PARENT);

     CONST_PTR(afw::geom::SkyWcs) wcs = exposure.getWcs();
     if (false) {
       //afw::geom::LinearTransform transform = wcs->linearizePixelToSky(center).getLinear();
       //PTR(afw::coord::Coord) skyCoord = wcs->pixelToSky(center);
       //afw::geom::Point2D sky(skyCoord->toIcrs().getRa(),skyCoord->toIcrs().getDec());
     }
     afw::geom::LinearTransform transform=afw::geom::LinearTransform::makeScaling(1.);

     double scaling = std::sqrt(std::abs(transform.computeDeterminant()));
     PTR(BFD::KWeight) kw;
     if (control.wIndex < 0) {
         kw.reset(new BFD::GaussianWeight(control.sigma*scaling));
     } else {
         kw.reset(new BFD::KSigmaWeight(control.sigma*scaling, control.wIndex));
     }


     Result result;
     try {
         result.flags[0] = false;
         double x0=center.getX();
         double y0=center.getY();

         PTR(BFD::PixelGalaxy<UseMoments>) pg;
         // Depending on how the noise is characterized pass different parameters to getPixelGalaxy
         if (control.useNoisePs) {
             PTR(daf::base::PropertyList) metadata = source.getTable()->getMetadata();
             kData = metadata->getArray<double>("kData");
             psData = metadata->getArray<double>("psData");
             pg = getPixelGalaxy(*kw, image, psf, transform, center, variance, x0, y0,
                                 control.reCentroid, control.reCentroidPsf, control.ignorePsf,
                                 control.interpOrder, kData, psData);
         } else if (control.useNoiseImagePs) {
             PTR(afw::image::Image<Pixel>) cov(new afw::image::Image<Pixel>(control.noiseImage));
             Table<> noisePs();
             std::vector<double> vec;
             PTR(BFD::PixelGalaxy<UseMoments>) pg = getPixelGalaxy(*kw, image, psf, transform, center, variance, x0, y0,
                                                                   control.reCentroid, control.reCentroidPsf, control.ignorePsf,
                                                                   control.interpOrder, vec, vec, cov);
         } else {
             pg = getPixelGalaxy(*kw, image, psf, transform, center, variance, x0, y0,
                                 control.reCentroid, control.reCentroidPsf, control.ignorePsf, control.interpOrder);
         }

         BFD::Moments<UseMoments> moments;
         BFD::MomentCovariance<UseMoments> momentsCov;

         DVector2 origin(0.);
         origin[0] = x0;
         origin[1] = y0;
         moments = pg->getMoments();
         momentsCov = pg->getCov();

         LOGL_DEBUG(trace3Logger,"Using center %g %g", origin[0], origin[1]);
         LOGL_DEBUG(trace4Logger,"got moments: %g", moments[0]);
         LOGL_DEBUG(trace4Logger,"             %g", moments[1]);
         LOGL_DEBUG(trace4Logger,"             %g", moments[2]);
         LOGL_DEBUG(trace4Logger,"             %g", moments[3]);
         LOGL_DEBUG(trace4Logger,"             %g", moments[4]);
         LOGL_DEBUG(trace4Logger,"             %g", moments[5]);
         LOGL_DEBUG(trace4Logger,"flux S/N: %g",moments[0]/std::sqrt(momentsCov(0,0)));
         LOGL_DEBUG(trace4Logger,"cov(0,0): %g",momentsCov(0,0));

         BFD::Moments<UseMoments> momentsPsf;
         if (!control.ignorePsf) {
             LOGL_DEBUG(trace3Logger,"Get moments of the psf");
             // Get the moments of the psf.  Use computeImage so that the psf will be centered on
             // the point x0,y0
             PTR(afw::image::Image<Pixel>) psfImage(
                 new afw::image::Image<Pixel>(*psf->computeImage(afw::geom::Point2D(x0,y0)), true));
             bool reCentroid = false;
             bool reCentroidPsf = false;
             bool ignorePsf = true;
             PTR(BFD::PixelGalaxy<UseMoments>)pgPsf = getPixelGalaxy(
                 *kw, psfImage, psf, transform, center, 1, x0, y0, reCentroid,
                 reCentroidPsf, ignorePsf, control.interpOrder);

             momentsPsf=pgPsf->getMoments();
             LOGL_DEBUG(trace4Logger,"Using center for PSF %g %g", x0, y0);
             LOGL_DEBUG(trace4Logger,"got moments for the PSF: %g", momentsPsf[0]);
             LOGL_DEBUG(trace4Logger,"                         %g", momentsPsf[1]);
             LOGL_DEBUG(trace4Logger,"                         %g", momentsPsf[2]);
             LOGL_DEBUG(trace4Logger,"                         %g", momentsPsf[3]);
             LOGL_DEBUG(trace4Logger,"                         %g", momentsPsf[4]);
             LOGL_DEBUG(trace4Logger,"                         %g", momentsPsf[5]);
         }

         if (control.shift) {
             double dx, dy;
             PTR(BFD::TemplateGalaxy<UseMoments>) shifted;
             try {
                 shifted.reset(BFD::newtonShift(*pg, dx, dy, control.shiftIterations));
             } catch(std::runtime_error const &error) {
                 LOGL_DEBUG(trace3Logger,"Shift failed");
                 result.flags[BfdKMomentResult::SHIFT_FAILED] = true;
                 throw;
             }

             if ( (dx*dx + dy*dy) > control.maxShift*control.maxShift) {
                 LOGL_DEBUG(trace3Logger,"Cenroid shift to null moments too large %g", std::sqrt(dx*dx + dy*dy) );
                 result.flags[BfdKMomentResult::SHIFT_LARGE] = true;
                 //result.flags[0] = true;
                 //return result;
             }
             LOGL_DEBUG(trace3Logger,"Shift to null centroid moments %g %g",dx,dy);
             origin[0] += dx;
             origin[1] += dy;
             moments = shifted->getMoments();
             momentsCov = shifted->getCov();

             LOGL_DEBUG(trace4Logger,"new moments: %g",moments[0]);
             LOGL_DEBUG(trace4Logger,"             %g",moments[1]);
             LOGL_DEBUG(trace4Logger,"             %g",moments[2]);
             LOGL_DEBUG(trace4Logger,"             %g",moments[3]);
             LOGL_DEBUG(trace4Logger,"             %g",moments[4]);
             LOGL_DEBUG(trace4Logger,"             %g",moments[5]);
             LOGL_DEBUG(trace4Logger,"flux S/N: %g",moments[0]/std::sqrt(momentsCov(0,0)));
             LOGL_DEBUG(trace4Logger,"cov(0,0): %g",momentsCov(0,0));

             if ( std::abs(moments[1]) > control.maxCentroid ||
                  std::abs(moments[2]) > control.maxCentroid ) {
                 LOGL_DEBUG(trace3Logger,"Shifted centroid moments too large %g %g, flagging",moments[1],moments[2]);
                 result.flags[BfdKMomentResult::SHIFT_CENTROID_LARGE] = true;
                 //result.flags[0] = true;
             }

             result.shift.setX(dx);
             result.shift.setY(dy);

         }
         result.center.setX(origin[0]);
         result.center.setY(origin[1]);

         if ( moments[0] < 0 ) {
             LOGL_DEBUG(trace3Logger,"Flux < 0, flagging");
             result.flags[BfdKMomentResult::FLUX_NEGATIVE] = true;
             result.flags[0] = true;
         }

         if ( moments[3] < 0 ) {
             LOGL_DEBUG(trace3Logger,"Size < 0, flagging");
             result.flags[BfdKMomentResult::SIZE_NEGATIVE] = true;
             result.flags[0] = true;
         }

         std::copy(moments.ptr(), moments.ptr()+NMoment, result.moments.getData());
         std::copy(momentsCov.ptr(), momentsCov.ptr()+NMoment*NMoment, result.momentsCov.getData());
         if (!control.ignorePsf) {
             std::copy(momentsPsf.ptr(), momentsPsf.ptr()+NMoment, result.momentsPsf.getData());
         }

     }  catch(std::runtime_error const &error) {
         LOGL_DEBUG(trace3Logger,"Runtime error: %s", error.what());
         result.flags[0] = true;
     } catch(...) {
         LOGL_DEBUG(trace3Logger,"Unknown exception, moment failed");
         result.flags[0] = true;
     }

     _resultKey.set(source, result);
     return;
 }



////////////////////////////////////////////////////////////////////////////////////////////
//-------------------- Prior Galaxy code --------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////
class PriorGalaxy::PriorGalaxyImpl
{
public:
    PriorGalaxyImpl(BfdKMomentControl const &ctrl):_ctrl(ctrl) {

        if (_ctrl.wIndex < 0) {
            kw.reset(new BFD::GaussianWeight(_ctrl.sigma));
        } else {
            kw.reset(new BFD::KSigmaWeight(_ctrl.sigma, _ctrl.wIndex));
        }

        //multiGal.reset(new BFD::MultiGalaxy<UseMoments>(*kw));
    }

    bool addImage(
             PTR(afw::image::Image<Pixel>)const & image,
             CONST_PTR(afw::detection::Psf) const & psf,
             afw::geom::LinearTransform const & transform,
             afw::geom::Point2D const & center,
             double noise) {
        LOG_LOGGER trace3Logger = LOG_GET("TRACE3.desc.bfd");
        //PTR(BFD::PixelGalaxy<UseMoments>) pg;
        try {
            double x0=center.getX();
            double y0=center.getY();
            LOGL_DEBUG(trace3Logger,"getting pixel galaxy for prior");
            pg = getPixelGalaxy(*kw, image, psf, transform, center, noise, x0, y0,
                                _ctrl.reCentroid, _ctrl.reCentroidPsf, _ctrl.ignorePsf, _ctrl.interpOrder);
            return true;
        } catch(const std::exception& ex) {

            std::cerr<<"Error adding prior template galaxy: "<<ex.what()<<std::endl;
            return false;
        }

        //if (pg) multiGal->add(pg.get());
        return false;
    }

    PTR(BFD::KWeight) kw;
    //PTR(BFD::MultiGalaxy<UseMoments>) multiGal;
    PTR(BFD::PixelGalaxy<UseMoments>) pg;
    BfdKMomentControl _ctrl;
};


PriorGalaxy::PriorGalaxy(BfdKMomentControl const &ctrl): _impl(new PriorGalaxyImpl(ctrl))
{
}

bool PriorGalaxy::addImage(
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure,
    afw::geom::Point2D const & center,
    double noise, bool useSrcFootprint)
{
    LOG_LOGGER trace3Logger = LOG_GET("TRACE3.desc.bfd");
    PTR(afw::detection::Psf const) psf = exposure.getPsf();
    PTR(afw::image::Image<Pixel>) imageExp = exposure.getMaskedImage().getImage();
    //PTR(afw::image::Wcs) wcs = exposure.getWcs();
    //afw::geom::LinearTransform transform = wcs->linearizePixelToSky(center).getLinear();
    afw::geom::LinearTransform transform;
    if (useSrcFootprint) {
        LOGL_DEBUG(trace3Logger,"Using footprint from src");
        PTR(afw::detection::Footprint) foot = source.getFootprint();
        afw::geom::Box2I box = foot->getBBox();
        PTR(afw::image::Image<Pixel>) image = std::make_shared<afw::image::Image<Pixel> >(*imageExp, box,
                                                                                            afw::image::PARENT);
        return _impl->addImage(image, psf, transform, center, noise);
    }
    return _impl->addImage(imageExp, psf, transform, center, noise);

}

bool PriorGalaxy::addScaledPSFImage(
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure,
    afw::geom::Point2D const & center,
    double noise, float scale)
{
    LOG_LOGGER trace3Logger = LOG_GET("TRACE3.desc.bfd");
    LOGL_DEBUG(trace3Logger,"entering addScaledPSFImage");
    PTR(afw::detection::Psf const) psf = exposure.getPsf();
    LOGL_DEBUG(trace3Logger,"creating psf image");
    PTR(afw::image::Image<Pixel>) psfImage(
                     new afw::image::Image<Pixel>(*psf->computeKernelImage(center), true));
    *psfImage *= scale;
    LOGL_DEBUG(trace3Logger,"scaling by %f",scale);
    afw::geom::LinearTransform transform;
    return _impl->addImage(psfImage, psf, transform, center, noise);

}

ndarray::Array<float, 1, 1> PriorGalaxy::getMoments() {
    ndarray::Array<float, 1, 1> result = ndarray::allocate(ndarray::makeVector(NMoment));
    //BFD::Moments<UseMoments> mom = _impl->multiGal->getMoments();
    BFD::Moments<UseMoments> mom = _impl->pg->getMoments();
    std::copy(mom.ptr(), mom.ptr() + NMoment, result.getData());
    // if (_impl->pixelGal) {
    //     BFD::Moments<UseMoments> mom = _impl->pixelGal->getMoments();
    //     std::copy(mom.ptr(), mom.ptr() + NMoment, result.getData());
    // }
    return result;
}

ndarray::Array<float,2,2>
PriorGalaxy::getMomentCov()
{
    //BFD::MomentCovariance<UseMoments> momCov = _impl->multiGal->getCov();
    BFD::MomentCovariance<UseMoments> momCov = _impl->pg->getCov();
    ndarray::Array<float, 2, 2> result = ndarray::allocate(ndarray::makeVector(NMoment,NMoment));
    std::copy(momCov.ptr(), momCov.ptr() + NMoment*NMoment, result.getData());
    return result;
}

ndarray::Array<float,2,2>
PriorGalaxy::getdMdG()
{
    BFD::MomentDerivs<UseMoments> deriv = _impl->pg->getDerivs();
    ndarray::Array<float, 2, 2> result = ndarray::allocate(ndarray::makeVector(NMoment,NMoment));
    std::copy(deriv.ptr(), deriv.ptr() + NMoment*NMoment, result.getData());
    return result;
}



////////////////////////////////////////////////////////////////////////////////////////////
//-------------------- Moment Prior code --------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////
class MomentPrior::MomentPriorImpl
{
public:
    MomentPriorImpl(double fluxMin, double fluxMax,
                    ndarray::Array<float,2,2> covariance,
                    bool invariantCovariance,
                    double noiseFactor,
                    double sigmaCutoff, double sigmaStep,
                    int nSamp, double sigmaBuffer,
                    bool selectOnly):ud(12454) {

        prior.reset(new BFD::KDTreePrior<UseMoments>(fluxMin, fluxMax, false, nSamp,
                                                     ud, sigmaBuffer, BFD::MID, false));

        BFD::MomentIndices<UseMoments>::MMatrix cov;
        std::copy(covariance.getData(), covariance.getData()+NMoment*NMoment, cov.ptr());
        BFD::MomentCovariance<UseMoments> momCov(cov);

        prior->setNominalCovariance(cov);
        prior->setInvariantCovariance(invariantCovariance);
        prior->setSamplingRange(sigmaCutoff, sigmaStep);
        //prior->setUniqueCeiling(200);
        if (selectOnly) prior->setSelectionOnly();
        if (noiseFactor > 1.) prior->addNoiseFactor(noiseFactor, ud);
    }

    PTR(BFD::KDTreePrior<UseMoments>) prior;
    ran::UniformDeviate ud;
}; 

MomentPrior::MomentPrior(double fluxMin, double fluxMax,
                         ndarray::Array<float,2,2> covariance,
                         bool invariantCovariance,
                         double noiseFactor,
                         double sigmaCutoff, double sigmaStep,
                         int nSamp, double sigmaBuffer, bool selectionOnly):
    _impl(new MomentPriorImpl(fluxMin, fluxMax, covariance, invariantCovariance, noiseFactor,
                              sigmaCutoff, sigmaStep, nSamp, sigmaBuffer, selectionOnly)),
    _fluxMin(fluxMin),
    _fluxMax(fluxMax),
    _init(true)
{}

void MomentPrior::addPriorGalaxy(
    PriorGalaxy const & gal, double maxXY, double weight, bool flip, long id)
{
    LOG_LOGGER trace3Logger = LOG_GET("TRACE3.desc.bfd");
    try {
        LOGL_DEBUG(trace3Logger,"adding template");

        double dx, dy;
        //PTR(BFD::TemplateGalaxy<UseMoments>) shifted(BFD::newtonShift(*gal._impl->multiGal, dx, dy,3));
        //PTR(BFD::TemplateGalaxy<UseMoments>) shifted(BFD::newtonShift(*gal._impl->pg, dx, dy,3));
        // LOGL_DEBUG(trace3Logger,"shifting template by %g %g",dx,dy);

        // debugLog.debug<2>("shifted moments: %g",shifted->getMoments()[0]);
        // debugLog.debug<2>("             %g",shifted->getMoments()[1]);
        // debugLog.debug<2>("             %g",shifted->getMoments()[2]);
        // debugLog.debug<2>("             %g",shifted->getMoments()[3]);
        // debugLog.debug<2>("             %g",shifted->getMoments()[4]);
        // debugLog.debug<2>("             %g",shifted->getMoments()[5]);

        if ( (dx*dx + dy*dy) > maxXY*maxXY) {
            LOGL_DEBUG(trace3Logger,"Cenroid shift to null moments too large %g, not adding galaxy", std::sqrt(dx*dx + dy*dy) );
            return;
        }
        //_impl->prior->addTemplate(*shifted, _impl->ud, maxXY, weight, flip, id);
        _impl->prior->addTemplate(*(gal._impl->pg), _impl->ud, maxXY, weight, flip, id);
    } catch(const std::exception& ex) {
        std::cerr<<"Error adding prior template galaxy: "<<ex.what()<<std::endl;
    }

}

void MomentPrior::prepare() {
    _impl->prior->prepare();
}


void MomentPrior::setMaxRatio(float ratio, int samples) {
    _impl->prior->setMaxRatio(ratio, samples);
}


BfdPqrResult
MomentPrior::getPqr(BfdKMomentResult const & moment) const
{
    BFD::MomentIndices<UseMoments>::MMatrix cov;
    std::copy(moment.momentsCov.getData(), moment.momentsCov.getData()+NMoment*NMoment, cov.ptr());
    BFD::MomentCovariance<UseMoments> momCov(cov);

    BFD::Moments<UseMoments> mom;
    std::copy(moment.moments.getData(), moment.moments.getData()+NMoment, mom.ptr());
    int nTemplates, nUnique;
    BFD::Pqr pqr = _impl->prior->getPqr2(mom, momCov, nTemplates, nUnique);

    BfdPqrResult result;
    std::copy(pqr.ptr(), pqr.ptr() + Npqr, result.pqr.getData());
    result.nTemplates = nTemplates;
    result.nUnique = nUnique;
    result.flags[0] = false;

    return result;
}

void
MomentPrior::getPqrCat(
    afw::table::SourceCatalog const & momentCat,
    afw::table::SourceCatalog &resultCat,
    int nthreads,
    int chunk
    ) const
{
    LOG_LOGGER trace3Logger = LOG_GET("TRACE3.desc.bfd");
#ifdef _OPENMP
    if (nthreads > 0) {
        omp_set_num_threads(nthreads);
    }
#endif

    BfdKMomentResultKey momentKey(momentCat.getSchema(),"bfd");
    BfdPqrKey pqrKey(resultCat.getSchema(),"bfd");
    afw::table::Key<afw::table::Flag> flagKey = resultCat.getSchema().find<afw::table::Flag>("bfd.flags").key;

 #ifdef _OPENMP
 #pragma omp parallel
     {
 #pragma omp for schedule(dynamic, chunk)
 #endif
         for(int i=0; i<momentCat.size(); ++i) {
             BfdKMomentResult moment;
#pragma omp critical
             {
                 moment = momentKey.get(momentCat[i]);
             }

              BfdPqrResult result;
              try {
                  result = this->getPqr(moment);
              }  catch(const std::exception& ex) {

#pragma omp critical
                  {
                      resultCat[i].set(flagKey, true);
                  }
              }

#pragma omp critical
              {
                  pqrKey.set(resultCat[i],result);
              }

         }
 #ifdef _OPENMP
     }
 #endif
}

BfdPqrResult
MomentPrior::getPqrArray(ndarray::Array<float,1,1> const & moments,
                    ndarray::Array<float,2,2> const & momentsCov) const
{
    BFD::MomentIndices<UseMoments>::MMatrix cov;
    std::copy(momentsCov.getData(), momentsCov.getData()+NMoment*NMoment, cov.ptr());
    BFD::MomentCovariance<UseMoments> momCov(cov);

    BFD::Moments<UseMoments> mom;
    std::copy(moments.getData(), moments.getData()+NMoment, mom.ptr());
    int nTemplates, nUnique;
    BFD::Pqr pqr = _impl->prior->getPqr2(mom, momCov, nTemplates, nUnique);

    // Take the negative log expansion here
    //BFD::Pqr nlPqr = pqr.neglog();
 
    BfdPqrResult result;
    std::copy(pqr.ptr(), pqr.ptr() + Npqr, result.pqr.getData());
    result.nTemplates = nTemplates;
    result.nUnique = nUnique;
    result.flags[0] = false;
    return result;
}


ndarray::Array<float,1,1>
MomentPrior::selectionProbability(ndarray::Array<float,2,2> const & covariance) const
{
     ndarray::Array<float,1,1> pqr = ndarray::allocate(ndarray::makeVector(NMoment));

     BFD::MomentIndices<UseMoments>::MMatrix cov;
     std::copy(covariance.getData(), covariance.getData()+NMoment*NMoment, cov.ptr());

     BFD::Pqr spqr = _impl->prior->selectionProbability(cov);
     BfdPqrResult result;
     std::copy(spqr.ptr(), spqr.ptr() + Npqr, pqr.getData());
     return pqr;
}

void
MomentPrior::writeFits(std::string name) {
    afw::table::SourceCatalog cat = this->getCatalog();
    cat.writeFits(name);
}

afw::table::SourceCatalog
MomentPrior::getCatalog() {

    LOG_LOGGER trace3Logger = LOG_GET("TRACE3.desc.bfd");
    afw::table::Schema schema = afw::table::SourceTable::makeMinimalSchema();
    BfdKMomentInfoResultKey key = BfdKMomentInfoResultKey::addFields(schema, "bfd");
    afw::table::SourceCatalog catalog(schema);

    const std::vector<BFD::MomentInfo<UseMoments> > & templates = _impl->prior->getTemplates();
    LOGL_DEBUG(trace3Logger,"prior has %d templates", templates.size());

    catalog.reserve(templates.size());

    for( int i=0; i<templates.size(); ++i) {
         BfdKMomentInfoResult result;
         std::copy(templates[i].m.cptr(),
                   templates[i].m.cptr()+NMoment,
                   result.moments.getData());
         std::copy(templates[i].deriv.cptr(),
                   templates[i].deriv.cptr()+NMoment*Npqr,
                   result.momentsDeriv.getData());
         result.weight = templates[i].weight;
         result.selectionWeight = templates[i].selectionWeight;
         result.id = templates[i].id;
         PTR(afw::table::SourceRecord) rec = catalog.addNew();
         key.set(*rec, result);
    }

    return catalog;
}

void
MomentPrior::readFits(std::string name, bool invarCov, float sample, int seed) {
    afw::table::SourceCatalog catalog = afw::table::SourceCatalog::readFits(name);
    this->addCatalog(catalog, invarCov);
}

double
MomentPrior::getTotalWeight() const {
    return _impl->prior->getTotalWeight();
}


void
MomentPrior::adjustTemplateWeights(std::vector<double> const &weights) {

    std::vector<BFD::MomentInfo<UseMoments> > & templates = _impl->prior->getTemplates();
    if (templates.size() != weights.size()) {
            throw LSST_EXCEPT(
                pex::exceptions::InvalidParameterError,
                "Cannot adjust weights with different number of templates");
    }

    for( int i=0; i<templates.size(); ++i) {
        templates[i].weight *= weights[i];
    }
}


std::vector<int>
MomentPrior::getTemplateWeights() const {

    std::vector<BFD::MomentInfo<UseMoments> > & templates = _impl->prior->getTemplates();
    std::vector<int> weights(templates.size());

    for( int i=0; i<templates.size(); ++i) weights[i] = templates[i].weight;

    return weights;
}

// copy BaseCatalog code
void
MomentPrior::addCatalog(afw::table::SourceCatalog const &catalog, bool invarCov, float sample, int seed) {

    PTR(daf::base::PropertyList) metadata = catalog.getTable()->getMetadata();
    float sampleWeight = 1./sample;

    if (!_init) {

        double fluxMin= metadata->getAsDouble("FLUXMIN");
        double fluxMax = metadata->getAsDouble("FLUXMAX");
        double varMin = metadata->getAsDouble("VARMIN");
        double varMax = metadata->getAsDouble("VARMAX");
        double totalWeight = metadata->getAsDouble("totalWeight");
        bool invariantCovariance = invarCov;//metadata->getAsBool("INVARIANTCOVARIANCE");

        std::vector<double> covVec = metadata->getArray<double>("COV");
        ndarray::Array<float,2,2> cov = ndarray::allocate(ndarray::makeVector(NMoment,NMoment));
        std::copy(covVec.begin(), covVec.begin() + NMoment*NMoment, cov.getData());

        double noiseFactor = metadata->getAsDouble("noiseFactor");
        double sigmaStep = metadata->getAsDouble("priorSigmaStep");
        int nSamp = metadata->getAsInt("NSAMPLE");
        double sigmaBuffer = metadata->getAsDouble("priorSigmaBuffer");
        double sigmaCutoff = metadata->getAsDouble("priorSigmaCutoff");
        bool selectionOnly = metadata->getAsBool("selectionOnly");
        _impl.reset(new MomentPriorImpl(fluxMin, fluxMax, cov, invariantCovariance,
                                        noiseFactor,sigmaCutoff, sigmaStep, nSamp,
                                        sigmaBuffer, selectionOnly));
        _impl->prior->setTotalWeight(totalWeight);
        _fluxMin = fluxMin;
        _fluxMax = fluxMax;
        _varMin = varMin;
        _varMax = varMax;
        _init = true;
    } else {
        double tol=0.1;
        // Need to do a more thourough check here that the prior can be added
        if (abs(_fluxMin-metadata->getAsDouble("FLUXMIN")) > tol ||
            abs(_fluxMax-metadata->getAsDouble("FLUXMAX")) > tol ||
            abs(_varMax-metadata->getAsDouble("VARMAX")) > tol ||
            abs(_varMin-metadata->getAsDouble("VARMIN")) > tol) {
            std::cout<<"BAD parameter"
                     <<abs(_fluxMin-metadata->getAsDouble("FLUXMIN"))<<" "
                     <<abs(_fluxMax-metadata->getAsDouble("FLUXMAX"))<<" "
                     <<abs(_varMin-metadata->getAsDouble("VARMIN"))<<" "
                     <<abs(_varMax-metadata->getAsDouble("VARMAX"))<<" "<<std::endl;
            throw LSST_EXCEPT(
                pex::exceptions::InvalidParameterError,
                "Cannot add priors with different parameters");
        }
        float currentWeight = _impl->prior->getTotalWeight();
        double catalogWeight = metadata->getAsDouble("totalWeight");
        _impl->prior->setTotalWeight(currentWeight+catalogWeight);

    }

    BfdKMomentInfoResultKey key(catalog.getSchema(),"bfd");
    ran::UniformDeviate myrand(seed);
    for (afw::table::SourceCatalog::const_iterator iter=catalog.begin(), end=catalog.end(); iter!=end; ++iter) {
        //if (rand > 0 && myrand() < sample) continue;
        if (myrand() > sample) continue;

        BfdKMomentInfoResult result = key.get(*iter);
        BFD::MomentInfo<UseMoments> info;
        info.weight = result.weight*sampleWeight;
        info.selectionWeight = result.selectionWeight*sampleWeight;
        info.id = result.id;
        std::copy(result.moments.getData(),
                  result.moments.getData() + NMoment,
                  info.m.ptr());
        std::copy(result.momentsDeriv.getData(),
                  result.momentsDeriv.getData() + NMoment*Npqr,
                  info.deriv.ptr());
        _impl->prior->getTemplates().push_back(info);
    }

}

void
MomentPrior::addCatalog(afw::table::BaseCatalog const &catalog, bool invarCov, float sample, int seed) {
 
    PTR(daf::base::PropertyList) metadata = catalog.getTable()->getMetadata();
    float sampleWeight = 1./sample;
    if (!_init) {

        double fluxMin= metadata->getAsDouble("FLUXMIN");
        double fluxMax = metadata->getAsDouble("FLUXMAX");
        double varMin = metadata->getAsDouble("VARMIN");
        double varMax = metadata->getAsDouble("VARMAX");
        double totalWeight = metadata->getAsDouble("totalWeight");
        bool invariantCovariance = invarCov;//metadata->getAsBool("INVARIANTCOVARIANCE");

        std::vector<double> covVec = metadata->getArray<double>("COV");
        ndarray::Array<float,2,2> cov = ndarray::allocate(ndarray::makeVector(NMoment,NMoment));
        std::copy(covVec.begin(), covVec.begin() + NMoment*NMoment, cov.getData());

        double noiseFactor = metadata->getAsDouble("noiseFactor");
        double sigmaStep = metadata->getAsDouble("priorSigmaStep");
        int nSamp = metadata->getAsInt("NSAMPLE");
        double sigmaBuffer = metadata->getAsDouble("priorSigmaBuffer");
        double sigmaCutoff = metadata->getAsDouble("priorSigmaCutoff");
        bool selectionOnly = metadata->getAsBool("selectionOnly");
        _impl.reset(new MomentPriorImpl(fluxMin, fluxMax, cov, invariantCovariance,
                                        noiseFactor,sigmaCutoff, sigmaStep, nSamp,
                                        sigmaBuffer, selectionOnly));
        _impl->prior->setTotalWeight(totalWeight);
        _fluxMin = fluxMin;
        _fluxMax = fluxMax;
        _varMin = varMin;
        _varMax = varMax;
        _init = true;
    } else {
        double tol=0.1;
        // Need to do a more thourough check here that the prior can be added
        if (abs(_fluxMin-metadata->getAsDouble("FLUXMIN")) > tol ||
            abs(_fluxMax-metadata->getAsDouble("FLUXMAX")) > tol ||
            abs(_varMax-metadata->getAsDouble("VARMAX")) > tol ||
            abs(_varMin-metadata->getAsDouble("VARMIN")) > tol) {
            std::cout<<"BAD parameter"
                     <<abs(_fluxMin-metadata->getAsDouble("FLUXMIN"))<<" "
                     <<abs(_fluxMax-metadata->getAsDouble("FLUXMAX"))<<" "
                     <<abs(_varMin-metadata->getAsDouble("VARMIN"))<<" "
                     <<abs(_varMax-metadata->getAsDouble("VARMAX"))<<" "<<std::endl;
            throw LSST_EXCEPT(
                pex::exceptions::InvalidParameterError,
                "Cannot add priors with different parameters");
        }
        float currentWeight = _impl->prior->getTotalWeight();
        double catalogWeight = metadata->getAsDouble("totalWeight");
        _impl->prior->setTotalWeight(currentWeight+catalogWeight);

    }

    BfdKMomentInfoResultKey key(catalog.getSchema(),"bfd");
    ran::UniformDeviate myrand(seed);
    for (afw::table::BaseCatalog::const_iterator iter=catalog.begin(), end=catalog.end(); iter!=end; ++iter) {
        //if (rand > 0 && myrand() < sample) continue;
        if (myrand() > sample) continue;

        BfdKMomentInfoResult result = key.get(*iter);
        BFD::MomentInfo<UseMoments> info;
        info.weight = result.weight*sampleWeight;
        info.selectionWeight = result.selectionWeight*sampleWeight;
        info.id = result.id;
        std::copy(result.moments.getData(),
                  result.moments.getData() + NMoment,
                  info.m.ptr());
        std::copy(result.momentsDeriv.getData(),
                  result.momentsDeriv.getData() + NMoment*Npqr,
                  info.deriv.ptr());
        _impl->prior->getTemplates().push_back(info);
    }

} 



////////////////////////////////////////////////////////////////////////////////////////////
//-------------------- Bfd KMoment Info Code --------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////
BfdKMomentInfoResultKey BfdKMomentInfoResultKey::addFields(
    afw::table::Schema & schema,
    std::string const & name
    )
{
    BfdKMomentInfoResultKey r;
    r._momentsKey = schema.addField<afw::table::Array<float> >(name+".moments",
                                                                    "k-space moments", NMoment);
    r._momentsDerivKey = schema.addField<afw::table::Array<float> >(name+".momentsDeriv",
                                                                    "k-space moments and derivatices", NMoment*Npqr);
    r._weightKey = schema.addField<float>(name+".weight","weight");
    r._selectionWeightKey = schema.addField<float>(name+".selectionWeight","selection weight");
    r._idKey = schema.addField<int>(name+".id","id");
    r._flagsKey[0] = schema.addField<afw::table::Flag>(name+".flags", "failure flag");
    r._flagsKey[1] = schema.addField<afw::table::Flag>(name+".flags.shift", "shift was too large");
    return r;

}

BfdKMomentInfoResultKey::BfdKMomentInfoResultKey(afw::table::Schema schema, std::string name) {
    _momentsKey = schema.find<afw::table::Array<float> >(name+".moments").key;
    _momentsDerivKey = schema.find<afw::table::Array<float> >(name+".momentsDeriv").key;
    _weightKey = schema.find<float>(name+".weight").key;
    _selectionWeightKey = schema.find<float>(name+".selectionWeight").key;
    _idKey = schema.find<int>(name+".id").key;
    _flagsKey[0] = schema.find<afw::table::Flag>(name+".flags").key;
    _flagsKey[1] = schema.find<afw::table::Flag>(name+".flags.shift").key;
}


void BfdKMomentInfoResultKey::set(afw::table::BaseRecord & record, BfdKMomentInfoResult const & value) const {
    record.set(_momentsKey, ndarray::flatten<1>(value.moments));
    record.set(_momentsDerivKey, ndarray::flatten<1>(value.momentsDeriv));
    record.set(_weightKey, value.weight);
    record.set(_selectionWeightKey, value.selectionWeight);
    record.set(_idKey, value.id);
    for (int b = 0; b < BfdKMomentInfoResult::N_FLAGS; ++b) {
          record.set(_flagsKey[b], value.flags[b]);
     }
}

BfdKMomentInfoResult BfdKMomentInfoResultKey::get(afw::table::BaseRecord const & record) const {
    BfdKMomentInfoResult result;
    ndarray::flatten<1>(result.moments) = record.get(_momentsKey);
    ndarray::flatten<1>(result.momentsDeriv) = record.get(_momentsDerivKey);
    result.weight = record.get(_weightKey);
    result.selectionWeight = record.get(_selectionWeightKey);
    result.id = record.get(_idKey);
    // Do not add flags for now because we need something better
         // for (int b = 0; b < BfdKMomentResult::N_FLAGS; ++b) {
         //   result.flags[b] = record.get(_flagsKey[b]);
         // }

    return result;
}



////////////////////////////////////////////////////////////////////////////////////////////
//-------------------- Bfd Pqr code --------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////
BfdPqrKey BfdPqrKey::addFields(
    afw::table::Schema & schema,
    std::string const & name
    )
{
    BfdPqrKey r;
    r._pqrKey = schema.addField<afw::table::Array<float> >(name+".pqr", "Pqr result", NMoment);
    r._nTemplatesKey = schema.addField<int>(name+".nTemplates", "number of templates");
    r._nUniqueKey = schema.addField<int>(name+".nUnique", "number of unique templates");
    r._flagsKey[0] = schema.addField<afw::table::Flag>(name+".flags", "failure flag");
    return r;
}

BfdPqrKey::BfdPqrKey(afw::table::Schema schema, std::string name) {
    _pqrKey = schema.find<afw::table::Array<float> >(name+".pqr").key;
    _nTemplatesKey = schema.find<int>(name+".nTemplates").key;
    _nUniqueKey = schema.find<int>(name+".nUnique").key;
    _flagsKey[0] = schema.find<afw::table::Flag>(name+".flags").key;
}


void BfdPqrKey::set(afw::table::BaseRecord & record, BfdPqrResult const & value) const {
    record.set(_pqrKey, value.pqr);
    record.set(_nTemplatesKey, value.nTemplates);
    record.set(_nUniqueKey, value.nUnique);
    for (int b = 0; b < BfdPqrResult::N_FLAGS; ++b) {
        record.set(_flagsKey[b], value.flags[b]);
    }
}

BfdPqrResult BfdPqrKey::get(afw::table::BaseRecord const & record) const {
    BfdPqrResult result;
    ndarray::flatten<1>(result.pqr) = record.get(_pqrKey);
    result.nTemplates = record.get(_nTemplatesKey);
    result.nUnique = record.get(_nUniqueKey);
    // Do not add flags for now because we need something better
    // for (int b = 0; b < BfdKMomentResult::N_FLAGS; ++b) {
    //   result.flags[b] = record.get(_flagsKey[b]);
    // }
    return result;
}
}}}

