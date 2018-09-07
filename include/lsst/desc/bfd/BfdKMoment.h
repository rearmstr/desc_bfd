
#ifndef LSST_DESC_BFD_BfdKMoment_h_INCLUDED
#define LSST_DESC_BFD_BfdMKoment_h_INCLUDED

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

#include "ndarray.h"

#include "lsst/pex/config.h"
#include "lsst/meas/base/exceptions.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/LinearTransform.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/desc/bfd/common.h"
#include "lsst/meas/base/Algorithm.h"

namespace lsst { namespace desc { namespace bfd {


class BfdKMomentResult;
class BfdKMomentControl;

class PriorGalaxy {
public:
    PriorGalaxy(BfdKMomentControl const &ctr);
    ~PriorGalaxy() = default;
    bool addImage(
        afw::table::SourceRecord & source,
        afw::image::Exposure<float> const & exposure,
        afw::geom::Point2D const & center,
        double noise, bool useSrcFootprint=true);

    bool addScaledPSFImage(
        afw::table::SourceRecord & source,
        afw::image::Exposure<float> const & exposure,
        afw::geom::Point2D const & center,
        double noise, float scale);


    ndarray::Array<float,1,1> getMoments();
    ndarray::Array<float,2,2> getMomentCov();
    ndarray::Array<float,2,2> getdMdG();

private:
    friend class MomentPrior;

    class PriorGalaxyImpl;
    PTR(PriorGalaxyImpl) _impl;
};

struct BfdPqrResult { 
     enum FlagBit { 
        FAILED=0, 
        N_FLAGS 
    }; 

    BfdPqrResult() { 
        pqr = ndarray::allocate(ndarray::makeVector(Npqr)); 
        flags[0] = true; 
        nTemplates = -1; 
        nUnique = -1; 
    } 

    ndarray::Array<float,1,1> pqr; 
    int nTemplates; 
    int nUnique; 

    std::bitset<N_FLAGS> flags; 
}; 


class BfdPqrKey { 
public: 
    static BfdPqrKey addFields( 
        afw::table::Schema & schema, 
        std::string const & name 
    ); 

    BfdPqrKey() {} 
    BfdPqrKey(afw::table::Schema schema, std::string name="bfd"); 

    virtual BfdPqrResult get(afw::table::BaseRecord const & record) const; 
    virtual void set(afw::table::BaseRecord & record, BfdPqrResult const & value) const; 

private: 

    afw::table::Key<afw::table::Array<float> > _pqrKey; 
    afw::table::Key<int> _nTemplatesKey; 
    afw::table::Key<int> _nUniqueKey; 
    afw::table::Key<afw::table::Flag> _flagsKey[BfdPqrResult::N_FLAGS]; 
}; 


class MomentPrior { 
public: 

    MomentPrior():_init(false) {}
    MomentPrior( 
        double fluxMin, double fluxMax, 
        ndarray::Array<float,2,2> covariance, 
        bool invariantCovariance, 
        double noiseFactor, 
        double sigmaCutoff, double sigmaStep, 
        int nSamp, double sigmaBuffer, 
        bool selectionOnly); 


    void addPriorGalaxy(PriorGalaxy const & gal, double maxXY, double weight, bool flip=false, long id=0); 

    BfdPqrResult getPqr(BfdKMomentResult const & moment) const; 
    void getPqrCat(afw::table::SourceCatalog const & momentCat, 
                   afw::table::SourceCatalog &resultCat, int nthreads=1, int chunk=100) const; 
    BfdPqrResult getPqrArray(ndarray::Array<float,1,1> const & moment, 
                                     ndarray::Array<float,2,2> const & momentCov) const; 

    ndarray::Array<float,1,1> selectionProbability(ndarray::Array<float,2,2> const & covariance) const; 

    void prepare(); 
    void writeFits(std::string name); 
    afw::table::SourceCatalog getCatalog(); 
    void readFits(std::string name, bool invarCov, float sample=-1., int seed=0); 
    void addCatalog(afw::table::BaseCatalog const &cat, bool invarCov, float sample=1., int seed=0); 
    void addCatalog(afw::table::SourceCatalog const &cat, bool invarCov, float sample=1., int seed=0); 

    double getFluxMin() const {return _fluxMin;} 
    double getFluxMax() const {return _fluxMax;} 
    void setVarianceLimits(double min, double max) { 
        _varMin = min; 
        _varMax = max; 
    } 
    double getVarMin() const {return _varMin;} 
    double getVarMax() const {return _varMax;} 

    double getTotalWeight() const; 

    void adjustTemplateWeights(std::vector<double> const &weights); 
    std::vector<int> getTemplateWeights() const; 

    void setMaxRatio(float ratio, int samples); 

 private: 
    class MomentPriorImpl; 
    PTR(MomentPriorImpl) _impl; 
    double _fluxMin; 
    double _fluxMax; 
    double _varMin; 
    double _varMax; 
    bool _init; 
}; 



struct BfdKMomentResult {

 enum FlagBit {
       FAILED=0,
       SHIFT_FAILED=1,
       SHIFT_LARGE=2,
       SHIFT_CENTROID_LARGE=3,
       TOO_BIG=4,
       VARIANCE=5,
       PARENT=6,
       FLUX_NEGATIVE=7,
       SIZE_NEGATIVE=8,
       SAT_CENTER=9,
       EMPTY_FOOTPRINT=10,
       N_FLAGS
   };

   BfdKMomentResult() {
       moments = ndarray::allocate(ndarray::makeVector(NMoment));
       momentsCov = ndarray::allocate(ndarray::makeVector(NMoment, NMoment));
       momentsPsf = ndarray::allocate(ndarray::makeVector(NMoment));
   }

   ndarray::Array<float,1,1> moments;
   ndarray::Array<float,1,1> momentsPsf;
   ndarray::Array<float,2,2> momentsCov;
   // Center used to measure the moments
   afw::geom::Point2D center;
   afw::geom::Point2D shift;
   std::bitset<N_FLAGS> flags;

};



class BfdKMomentResultKey {
public:
   static BfdKMomentResultKey addFields(
       afw::table::Schema & schema,
       std::string const & name
   );

   BfdKMomentResultKey() {}
   BfdKMomentResultKey(afw::table::Schema schema, std::string name="bfd");

   virtual BfdKMomentResult get(afw::table::BaseRecord const & record) const;
   virtual void set(afw::table::BaseRecord & record, BfdKMomentResult const & value) const;


private:

   afw::table::Key<afw::table::Array<float> > _momentsKey;
   afw::table::Key<afw::table::Array<float> > _momentsCovKey;
   afw::table::Key<afw::table::Array<float> > _momentsPsfKey;
   afw::table::PointKey<double> _centerKey;
   afw::table::PointKey<double> _shiftKey;
   afw::table::Key<afw::table::Flag> _flagsKey[BfdKMomentResult::N_FLAGS];


};

struct BfdKMomentControl {

   LSST_CONTROL_FIELD(
       sigma, double,
       "Sigma of weight funciton"
       );

   LSST_CONTROL_FIELD(
       wIndex, double,
       "Power of radial weight in (1-(k*sigma)**2)**n. If < 0 then the weight "
       "function is a Gaussian");

   LSST_CONTROL_FIELD(
       ignorePsf, bool,
       "Ignore PSF when calculating moments and just use a delta function. "
       "This is useful for debugging purposes."
       );

   LSST_CONTROL_FIELD(
       shift, bool,
       "Shift moments so that the centroid moments are zero."
       );

   LSST_CONTROL_FIELD(
       reCentroid, bool,
       "Find centroid of image, ignoring centroid information given."
       );

   LSST_CONTROL_FIELD(
       reCentroidPsf, bool,
       "Find centroid of psf, and shift to center"
       );

   LSST_CONTROL_FIELD(
       shiftIterations, int,
       "Number gauss-newton iterations to find zero centroid moments");

   LSST_CONTROL_FIELD(
       maxShift, double,
       "maximum shift allowed without flagging as problem");

   LSST_CONTROL_FIELD(
       maxCentroid, double,
       "maximum centroid allowed without flagging");

   LSST_CONTROL_FIELD(
       calculateVariance, bool,
       "if true it will calculate the noise from the image");

   LSST_CONTROL_FIELD(
       useTableVariance, bool,
       "if true it will get the noise from the metadata");

   LSST_CONTROL_FIELD(
       useRecVariance, bool,
       "if true it will get the noise from the SourceRecord");

   LSST_CONTROL_FIELD(
       interpOrder, int,
       "PSF Interpolation order");

   LSST_CONTROL_FIELD(
       useNoisePs, bool,
       "use noise power spectrum from the metadata");

    LSST_CONTROL_FIELD(
       useNoiseImagePs, bool,
       "use noise power spectrum from an image");

   LSST_CONTROL_FIELD(
       noiseImage, std::string,
       "location of noise image");

   BfdKMomentControl() :
       sigma(1.),
       wIndex(-1.),
       ignorePsf(false),
       shift(false),
       reCentroid(true),
       reCentroidPsf(false),
       shiftIterations(2),
       maxShift(3),
       maxCentroid(0.01),
       calculateVariance(false),
       useTableVariance(false),
       useRecVariance(true),
       interpOrder(5),
       useNoisePs(false),
       useNoiseImagePs(false),
       noiseImage("")
   {}
};



 class BfdKMoment : public meas::base::SimpleAlgorithm  {
public:

   typedef BfdKMomentControl Control;
   typedef BfdKMomentResult Result;
   typedef BfdKMomentResultKey ResultKey;

   BfdKMoment(
	      BfdKMomentControl const & ctrl,
	       std::string const & name,
	       afw::table::Schema & schema
	       );

   static Result measure_image(
        BfdKMomentControl const & ctrl,
        PTR(afw::image::Image<Pixel>) const & image,
        CONST_PTR(afw::detection::Psf) const & psf,
        afw::geom::LinearTransform const & transform,
        afw::geom::Point2D const & center,
        double noise=1.
    );

   static Result measure_exp(
       BfdKMomentControl const & ctrl,
       afw::geom::Box2I box,
       afw::image::Exposure<Pixel> const & exposure,
       afw::geom::Point2D const & center,
       double variance
    );

    void measure(
       afw::table::SourceRecord & measRecord,
       afw::image::Exposure<Pixel> const & exposure
     ) const;

    void fail(
	 afw::table::SourceRecord & measRecord,
	 meas::base::MeasurementError * error
     ) const;
   
private:

    //friend class BfdMomentControl;

     void _apply(
       afw::table::SourceRecord & source,
       afw::image::Exposure<Pixel> const & exposure,
       afw::geom::Point2D const & center
       ) const;

   Control _ctrl;
   ResultKey _resultKey;
   mutable std::map<int,float> _noiseCache;
};

struct BfdKMomentInfoResult { 
     enum FlagBit { 
        FAILED=0, 
        N_FLAGS 
    }; 

    BfdKMomentInfoResult() { 
        moments = ndarray::allocate(ndarray::makeVector(NMoment)); 
        momentsDeriv = ndarray::allocate(ndarray::makeVector(NMoment, Npqr)); 
    } 

    ndarray::Array<float,1,1> moments; 
    ndarray::Array<float,2,2> momentsDeriv; 
    int id;// make this larger - correspond to detection id? 
    float weight; 
    float selectionWeight; 

    std::bitset<N_FLAGS> flags; 

}; 


class BfdKMomentInfoResultKey {
public:
    static BfdKMomentInfoResultKey addFields(
        afw::table::Schema & schema,
        std::string const & name
    );

    BfdKMomentInfoResultKey() {}
    BfdKMomentInfoResultKey(afw::table::Schema schema, std::string name="bfd");

    virtual BfdKMomentInfoResult get(afw::table::BaseRecord const & record) const;
    virtual void set(afw::table::BaseRecord & record, BfdKMomentInfoResult const & value) const;


private: 

    afw::table::Key<afw::table::Array<float> > _momentsKey; 
    afw::table::Key<afw::table::Array<float> > _momentsDerivKey; 
    afw::table::Key<float> _weightKey; 
    afw::table::Key<float> _selectionWeightKey; 
    afw::table::Key<int> _idKey; 
    afw::table::Key<afw::table::Flag> _flagsKey[BfdKMomentInfoResult::N_FLAGS]; 


}; 

//
//
//struct BfdKMomentInfoControl : public algorithms::AlgorithmControl { 
//
//    LSST_CONTROL_FIELD( 
//        sigma, double, 
//        "Sigma of weight funciton" 
//        ); 
//
//    LSST_CONTROL_FIELD( 
//        wIndex, double, 
//        "Power of radial weight in (1-(k*sigma)**2)**n. If < 0 then the weight " 
//        "function is a Gaussian"); 
//
//    LSST_CONTROL_FIELD( 
//        ignorePsf, bool, 
//        "Ignore PSF when calculating moments and just use a delta function. " 
//        "This is useful for debugging purposes." 
//        ); 
//
//    LSST_CONTROL_FIELD( 
//        shift, bool, 
//        "Shift moments so that the centroid moments are zero." 
//        ); 
//
//    LSST_CONTROL_FIELD( 
//        reCentroid, bool, 
//        "Find centroid of image, ignoring centroid information given." 
//        ); 
//
//    LSST_CONTROL_FIELD( 
//        reCentroidPsf, bool, 
//        "Find centroid of psf and shift." 
//        ); 
//
//    LSST_CONTROL_FIELD( 
//        shiftIterations, int, 
//        "Number gauss-newton iterations to find zero centroid moments"); 
//
//    LSST_CONTROL_FIELD( 
//        maxShift, double, 
//        "maximum shift allowed without flagging as problem"); 
//
//    LSST_CONTROL_FIELD( 
//        fluxMax, float, 
//        "maximum flux considered"); 
//
//    LSST_CONTROL_FIELD( 
//        fluxMin, float, 
//        "minimum flux considered"); 
//
//    LSST_CONTROL_FIELD( 
//        sigmaStep, float, 
//        "step in sigma for generating templates"); 
//
//    LSST_CONTROL_FIELD( 
//        sigmaCutoff, float, 
//        "maximum sigma used for keeping templates from nominal covariance matrix"); 
//
//    LSST_CONTROL_FIELD( 
//        noiseFactor, float, 
//        "add Noise to moments"); 
//
//    LSST_CONTROL_FIELD( 
//        doDeselection, bool, 
//        "only do deselection"); 
//
//    LSST_CONTROL_FIELD( 
//        interpOrder, int, 
//        "PSF Interpolation order"); 
//
//    // Not sure if array can be field 
//    // LSST_CONTROL_FIELD( 
//    //     nominalCovariance, ndarray::Array<float,2,2>, 
//    //     "nominal covariance used in cutoff"); 
//    ndarray::Array<float,2,2> nominalCov; 
//
//    BfdKMomentInfoControl() : 
//        algorithms::AlgorithmControl("bfd.kmomentsInfo", 2.5), 
//        sigma(1.), 
//        wIndex(-1.), 
//        ignorePsf(false), 
//        shift(false), 
//        reCentroid(false), 
//        reCentroidPsf(false), 
//        shiftIterations(2), 
//        maxShift(3), 
//        fluxMin(5.), 
//        fluxMax(1000.), 
//        sigmaStep(1.), 
//        sigmaCutoff(5.5), 
//        noiseFactor(1.), 
//        doDeselection(false), 
//        interpOrder(5) 
//        { 
//            nominalCov = ndarray::allocate(ndarray::makeVector(NMoment,NMoment)); 
//            for(int i=0; i<NMoment; ++i) nominalCov[i][i]=1; 
//        } 
//private: 
//    virtual PTR(algorithms::AlgorithmControl) _clone() const; 
//    virtual PTR(algorithms::Algorithm) _makeAlgorithm( 
//        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata 
//    ) const; 
//}; 
//
//
//
//
//class BfdKMomentInfo : public algorithms::Algorithm { 
//public: 
//
//    typedef BfdKMomentInfoControl Control; 
//    typedef BfdKMomentInfoResult Result; 
//    typedef BfdKMomentInfoResultKey ResultKey; 
//
//    BfdKMomentInfo(BfdKMomentInfoControl const & ctrl, afw::table::Schema & schema); 
//
//    void writeFits(std::string name) { _infoCatalog.writeFits(name);} 
//
//     // static Result measure( 
//     //     BfdKMomentInfoControl const & ctrl, 
//     //     PTR(afw::image::Image<Pixel>) const & image, 
//     //     CONST_PTR(afw::detection::Psf) const & psf, 
//     //     afw::geom::LinearTransform const & transform, 
//     //     afw::geom::Point2D const & center, 
//     //     double noise=1. 
//     //    ); 
//private: 
//
//    template <typename PixelT> 
//    void _apply( 
//        afw::table::SourceRecord & source, 
//        afw::image::Exposure<PixelT> const & exposure, 
//        afw::geom::Point2D const & center 
//        ) const; 
//
//    Control _ctrl; 
//    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(BfdKMomentInfo); 
//    ResultKey _resultKey; 
//    afw::table::Schema _schema; 
//    mutable afw::table::SourceCatalog _infoCatalog; 
//
//    class PriorImpl; 
//    PTR(PriorImpl) _impl; 
//}; 
//
//
//
//struct BfdPixelGalaxyResult { 
//     enum FlagBit { 
//        FAILED=0, 
//        N_FLAGS 
//    }; 
//
//    BfdPixelGalaxyResult(int _size): 
//    size(_size), 
//    d2k(0.), 
//    dcIndex(-1), 
//    noise(1.) 
//        { 
//            reGal = ndarray::allocate(ndarray::makeVector(size)); 
//            imGal = ndarray::allocate(ndarray::makeVector(size)); 
//            rePsf = ndarray::allocate(ndarray::makeVector(size)); 
//            imPsf = ndarray::allocate(ndarray::makeVector(size)); 
//            kx = ndarray::allocate(ndarray::makeVector(size)); 
//            ky = ndarray::allocate(ndarray::makeVector(size)); 
//        } 
//
//    int size; 
//    ndarray::Array<float,1,1> reGal; 
//    ndarray::Array<float,1,1> imGal; 
//    ndarray::Array<float,1,1> rePsf; 
//    ndarray::Array<float,1,1> imPsf; 
//    ndarray::Array<float,1,1> kx; 
//    ndarray::Array<float,1,1> ky; 
//    float d2k; // area of grid cell in k-space 
//    int dcIndex; // vector index of point at zero frequency 
//    float noise; 
//
//    // Center used to measure the moments 
//    afw::geom::Point2D center; 
//#ifndef SWIG 
//    std::bitset<N_FLAGS> flags; 
//#endif 
//
//}; 
//
//
//class BfdPixelGalaxyResultKey { 
//public: 
//    static BfdPixelGalaxyResultKey addFields( 
//        afw::table::Schema & schema, 
//        std::string const & name 
//    ); 
//
//    BfdPixelGalaxyResultKey() {} 
//    BfdPixelGalaxyResultKey(afw::table::Schema schema, std::string name="bfd"); 
//    BfdPixelGalaxyResultKey( 
//        afw::table::Key<int> const & sizeKey, 
//        afw::table::Key<afw::table::Array<float> > const & reGalKey, 
//        afw::table::Key<afw::table::Array<float> > const & imGalKey, 
//        afw::table::Key<afw::table::Array<float> > const & rePsfKey, 
//        afw::table::Key<afw::table::Array<float> > const & imPsfKey, 
//        afw::table::Key<afw::table::Array<float> > const & kxKey, 
//        afw::table::Key<afw::table::Array<float> > const & kyKey, 
//        afw::table::Key<float> const & d2kKey, 
//        afw::table::Key<int> const & dcIndexKey, 
//        afw::table::Key<float> const & noiseKey, 
//        afw::table::Key<afw::table::Point<double> > const & centerKey) 
//        : _sizeKey(sizeKey), 
//          _reGalKey(reGalKey), 
//          _imGalKey(imGalKey), 
//          _rePsfKey(rePsfKey), 
//          _imPsfKey(imPsfKey), 
//          _kxKey(kxKey), 
//          _kyKey(kyKey), 
//          _d2kKey(d2kKey), 
//          _dcIndexKey(dcIndexKey), 
//          _noiseKey(noiseKey), 
//          _centerKey(centerKey) {} 
//
//
//
//    virtual BfdPixelGalaxyResult get(afw::table::BaseRecord const & record) const; 
//    virtual void set(afw::table::BaseRecord & record, BfdPixelGalaxyResult const & value) const; 
//
//    //ndarray1::Array<float,2,2> getCov(); 
//private: 
//
//    afw::table::Key<int> _sizeKey; 
//    afw::table::Key<afw::table::Array<float> > _reGalKey; 
//    afw::table::Key<afw::table::Array<float> > _imGalKey; 
//    afw::table::Key<afw::table::Array<float> > _rePsfKey; 
//    afw::table::Key<afw::table::Array<float> > _imPsfKey; 
//    afw::table::Key<afw::table::Array<float> > _kxKey; 
//    afw::table::Key<afw::table::Array<float> > _kyKey; 
//    afw::table::Key<float> _d2kKey; 
//    afw::table::Key<int> _dcIndexKey; 
//    afw::table::Key<float> _noiseKey; 
//    afw::table::Key<afw::table::Point<double> > _centerKey; 
//    afw::table::Key<afw::table::Flag> _flagsKey[BfdPixelGalaxyResult::N_FLAGS]; 
//
//}; 
//
//
//
//struct BfdPixelGalaxyControl : public algorithms::AlgorithmControl { 
//
//    LSST_CONTROL_FIELD( 
//        sigma, double, 
//        "Sigma of weight funciton" 
//        ); 
//
//    LSST_CONTROL_FIELD( 
//        wIndex, double, 
//        "Power of radial weight in (1-(k*sigma)**2)**n. If < 0 then the weight " 
//        "function is a Gaussian"); 
//
//    BfdPixelGalaxyControl() : 
//        algorithms::AlgorithmControl("bfd.PixelGalaxy", 2.5), 
//        sigma(1.), 
//        wIndex(-1.) 
//    {} 
//
//private: 
//    virtual PTR(algorithms::AlgorithmControl) _clone() const; 
//    virtual PTR(algorithms::Algorithm) _makeAlgorithm( 
//        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata 
//    ) const; 
//}; 
//
//
//class BfdPixelGalaxy : public algorithms::Algorithm { 
//public: 
//
//    typedef BfdPixelGalaxyControl Control; 
//    typedef BfdPixelGalaxyResult Result; 
//    typedef BfdPixelGalaxyResultKey ResultKey; 
//
//    BfdPixelGalaxy(BfdPixelGalaxyControl const & ctrl, afw::table::Schema & schema); 
//
//private: 
//
//    template <typename PixelT> 
//    void _apply( 
//        afw::table::SourceRecord & source, 
//        afw::image::Exposure<PixelT> const & exposure, 
//        afw::geom::Point2D const & center 
//        ) const; 
//
//    Control _ctrl; 
//    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(BfdPixelGalaxy); 
//    ResultKey _resultKey; 
//    mutable std::map<int,float> _noiseCache; 
//}; 
//
//
//
}}}
#endif
//
