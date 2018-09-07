// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#include "pybind11/pybind11.h"
#include "ndarray/pybind11.h"

#include "lsst/pex/config/python.h"
#include "lsst/desc/bfd/BfdKMoment.h"
#include "lsst/afw/detection/Psf.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace desc {
namespace bfd {

// Custom wrapper for views to std::bitset.
 template <int N>
 class BitSetView {
 public:
     explicit BitSetView(std::bitset<N> const *target) : _target(target) {}

     bool operator[](int i) const { return (*_target)[i]; }

     constexpr std::size_t size() const { return N; }

     template <typename PyParent>
     static void declare(PyParent &parent) {
         py::class_<BitSetView<N>> cls(parent, ("BitSetView" + std::to_string(N)).c_str());
         cls.def("__getitem__", &BitSetView<N>::operator[]);
         cls.def("__len__", &BitSetView<N>::size);
     }

 private:
     std::bitset<N> const *_target;
 };

PYBIND11_MODULE(bfdkmoment, mod) {

     py::class_<BfdKMoment, std::shared_ptr<BfdKMoment>, meas::base::SimpleAlgorithm> clsBfdKMoment(mod, "BfdKMoment");
     clsBfdKMoment.def(py::init<BfdKMomentControl const &, std::string const &,
                                afw::table::Schema &>(),
                       "ctrl"_a, "name"_a, "schema"_a);
    clsBfdKMoment.def_static("measure_image", &BfdKMoment::measure_image, "ctrl"_a, "image"_a, "psf"_a, "transform"_a, "center"_a, "noise"_a);
    clsBfdKMoment.def_static("measure_exp", &BfdKMoment::measure_exp, "ctrl"_a, "box"_a, "exposure"_a, "center"_a, "variance"_a);
    clsBfdKMoment.def("measure", &BfdKMoment::measure, "measRecord"_a, "exposure"_a);
    clsBfdKMoment.def("fail", &BfdKMoment::fail, "measRecord"_a, "error"_a);

    py::class_<BfdKMomentResult, std::shared_ptr<BfdKMomentResult>> clsBfdKMomentResult(mod, "BfdKMomentResult");
    clsBfdKMomentResult.def(py::init<>());
    clsBfdKMomentResult.attr("FAILED") = py::cast(int(BfdKMomentResult::FAILED));
    clsBfdKMomentResult.attr("SHIFT_FAILED") = py::cast(int(BfdKMomentResult::SHIFT_FAILED));
    clsBfdKMomentResult.attr("SHIFT_LARGE") = py::cast(int(BfdKMomentResult::SHIFT_LARGE));
    clsBfdKMomentResult.attr("SHIFT_CENTROID_LARGE") = py::cast(int(BfdKMomentResult::SHIFT_CENTROID_LARGE));
    clsBfdKMomentResult.attr("TOO_BIG") = py::cast(int(BfdKMomentResult::TOO_BIG));
    clsBfdKMomentResult.attr("VARIANCE") = py::cast(int(BfdKMomentResult::VARIANCE));
    clsBfdKMomentResult.attr("PARENT") = py::cast(int(BfdKMomentResult::PARENT));
    clsBfdKMomentResult.attr("FLUX_NEGATIVE") = py::cast(int(BfdKMomentResult::FLUX_NEGATIVE));
    clsBfdKMomentResult.attr("SIZE_NEGATIVE") = py::cast(int(BfdKMomentResult::SIZE_NEGATIVE));
    clsBfdKMomentResult.attr("SAT_CENTER") = py::cast(int(BfdKMomentResult::SAT_CENTER));
    clsBfdKMomentResult.attr("EMPTY_FOOTPRINT") = py::cast(int(BfdKMomentResult::EMPTY_FOOTPRINT));

    clsBfdKMomentResult.def_readonly("moments", &BfdKMomentResult::moments);
    clsBfdKMomentResult.def_readonly("momentsCov", &BfdKMomentResult::momentsCov);
    clsBfdKMomentResult.def_readonly("momentsPsf", &BfdKMomentResult::momentsPsf);
    clsBfdKMomentResult.def_readonly("center", &BfdKMomentResult::center);
    clsBfdKMomentResult.def_readonly("shift", &BfdKMomentResult::shift);

    clsBfdKMomentResult.def_property_readonly(
        		      "flags", [](BfdKMomentResult const &self) { return BitSetView<BfdKMomentResult::N_FLAGS>(&self.flags); },
        		      py::return_value_policy::reference_internal);

    py::class_<BfdKMomentControl, std::shared_ptr<BfdKMomentControl>> clsBfdKMomentControl(mod, "BfdKMomentControl");
    clsBfdKMomentControl.def(py::init<>());
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, sigma);
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, wIndex);
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, ignorePsf);
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, shift);
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, reCentroid);
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, reCentroidPsf);
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, shiftIterations);
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, maxShift);
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, maxCentroid);
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, calculateVariance);
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, useTableVariance);
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, useRecVariance);
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, interpOrder);
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, useNoisePs);
    LSST_DECLARE_CONTROL_FIELD(clsBfdKMomentControl, BfdKMomentControl, noiseImage);


    py::class_<BfdPqrResult, std::shared_ptr<BfdPqrResult>> clsBfdPqrResult(mod, "BfdPqrResult");
    clsBfdPqrResult.def(py::init<>());
    clsBfdPqrResult.attr("FAILED") = py::cast(int(BfdPqrResult::FAILED));

    clsBfdPqrResult.def_readonly("pqr", &BfdPqrResult::pqr);
    clsBfdPqrResult.def_readonly("nTemplates", &BfdPqrResult::nTemplates);
    clsBfdPqrResult.def_readonly("nUnique", &BfdPqrResult::nUnique);

    py::class_<MomentPrior, std::shared_ptr<MomentPrior>> clsMomentPrior(mod, "MomentPrior");
    clsMomentPrior.def(py::init<>());
    clsMomentPrior.def(py::init<double, double, ndarray::Array<float,2,2>, bool, double, double, double,
                       int, double, bool>());
    clsMomentPrior.def("addPriorGalaxy", &MomentPrior::addPriorGalaxy, "gal"_a, "maxXY"_a, "weight"_a, "flip"_a, "id"_a);
    clsMomentPrior.def("getPqr", &MomentPrior::getPqr, "moment"_a);
    clsMomentPrior.def("getPqrCat", &MomentPrior::getPqrCat, "momentCat"_a, "resultCat"_a, "nthreads"_a, "chunk"_a);
    clsMomentPrior.def("getPqrArray", &MomentPrior::getPqrArray, "moment"_a, "momentCov"_a);
    clsMomentPrior.def("selectionProbability", &MomentPrior::selectionProbability, "covariance"_a);
    clsMomentPrior.def("prepare", &MomentPrior::prepare);
    clsMomentPrior.def("writeFits", &MomentPrior::writeFits, "name"_a);
    clsMomentPrior.def("getCatalog", &MomentPrior::getCatalog);
    clsMomentPrior.def("addCatalog",
                        py::overload_cast<afw::table::SourceCatalog const &, bool, float, int>(&MomentPrior::addCatalog),
                        "cat"_a, "invarCov"_a, "sample"_a, "seed"_a
                        );
    clsMomentPrior.def("addCatalog",
                        py::overload_cast<afw::table::BaseCatalog const &, bool, float, int>(&MomentPrior::addCatalog),
                        "cat"_a, "invarCov"_a, "sample"_a, "seed"_a
                        );
    clsMomentPrior.def("getFluxMin", &MomentPrior::getFluxMin);
    clsMomentPrior.def("getFluxMax", &MomentPrior::getFluxMax);
    clsMomentPrior.def("getVarMin", &MomentPrior::getVarMin);
    clsMomentPrior.def("getVarMax", &MomentPrior::getVarMax);
    clsMomentPrior.def("setVarianceLimits", &MomentPrior::setVarianceLimits, "min"_a, "max"_a);
    clsMomentPrior.def("getTotalWeight", &MomentPrior::getTotalWeight);
    clsMomentPrior.def("adjustTemplateWeights", &MomentPrior::adjustTemplateWeights, "weights"_a);
    clsMomentPrior.def("getTemplateWeights", &MomentPrior::getTemplateWeights);
    clsMomentPrior.def("setMaxRatio", &MomentPrior::setMaxRatio, "ratio"_a, "samples"_a);


    py::class_<BfdKMomentInfoResult, std::shared_ptr<BfdKMomentInfoResult>> clsBfdKMomentInfoResult(mod, "BfdKMomentInfoResult");
    clsBfdKMomentInfoResult.def(py::init<>());
    clsBfdKMomentInfoResult.attr("FAILED") = py::cast(int(BfdKMomentInfoResult::FAILED));

    clsBfdKMomentInfoResult.def_readonly("moments", &BfdKMomentInfoResult::moments);
    clsBfdKMomentInfoResult.def_readonly("momentsDeriv", &BfdKMomentInfoResult::momentsDeriv);
    clsBfdKMomentInfoResult.def_readonly("weight", &BfdKMomentInfoResult::weight);
    clsBfdKMomentInfoResult.def_readonly("selectionWeight", &BfdKMomentInfoResult::selectionWeight);

    py::class_<PriorGalaxy, std::shared_ptr<PriorGalaxy>> clsPriorGalaxy(mod, "PriorGalaxy");
    clsPriorGalaxy.def(py::init<BfdKMomentControl const &>(), "ctr"_a);
    clsPriorGalaxy.def("addImage",&PriorGalaxy::addImage, "source"_a, "exposure"_a, "center"_a, "noise"_a, "scale"_a);
    clsPriorGalaxy.def("addScaledPSFImage",&PriorGalaxy::addScaledPSFImage);
    clsPriorGalaxy.def("getMoments",&PriorGalaxy::getMoments);
    clsPriorGalaxy.def("getMomentCov",&PriorGalaxy::getMomentCov);
    clsPriorGalaxy.def("getdMdG",&PriorGalaxy::getdMdG);

    py::class_<BfdPqrKey, std::shared_ptr<BfdPqrKey>> clsBfdPqrKey(mod, "BfdPqrKey");
    clsBfdPqrKey.def(py::init<>());
    clsBfdPqrKey.def(py::init<afw::table::Schema, std::string>(), "schema"_a, "name"_a);
    clsBfdPqrKey.def_static("addFields", &BfdPqrKey::addFields, "schema"_a, "name"_a);
    clsBfdPqrKey.def("get", &BfdPqrKey::get, "record"_a);
    clsBfdPqrKey.def("set", &BfdPqrKey::set, "record"_a, "value"_a);
    


}

}
}
}
