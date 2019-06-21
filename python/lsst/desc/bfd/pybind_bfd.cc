#include "pybind11/eigen.h"
#include "pybind11/operators.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "BfdConfig.h"
#include "Galaxy.h"
#include "KGalaxy.h"
#include "KdPrior.h"
#include "Moment.h"
#include "Pqr.h"
#include "Random.h"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace bfd;

namespace lsst {
namespace desc {
namespace bfd {

namespace {

std::string optionsToString(std::string name, bool cen, bool conc, bool mag,
                            int colors, bool ufloat) {
  std::string label(name);
  name += "_";
  if (cen)
    label += "1";
  else
    label += "0";

  if (conc)
    label += "1";
  else
    label += "0";

  if (mag)
    label += "1";
  else
    label += "0";

  label += std::to_string(colors);

  if (ufloat)
    label += "F";
  else
    label += "D";

  return label;
}

template <bool FIX_CENTER, bool USE_CONC, bool USE_MAG, int N_COLORS,
          bool USE_FLOAT>
using PyBfdConfig =
    py::class_<BfdConfig<FIX_CENTER, USE_CONC, USE_MAG, N_COLORS, USE_FLOAT>>;

template <class CONFIG>
using PyMoment = py::class_<Moment<CONFIG>>;

template <class CONFIG>
using PyMomentCov = py::class_<MomentCovariance<CONFIG>>;

template <class FP>
using PyKSigmaWeight = py::class_<KSigmaWeight<FP>>;

template <class CONFIG>
using PyKGalaxy = py::class_<KGalaxy<CONFIG>>;

template <class CONFIG>
using PyTargetGalaxy = py::class_<TargetGalaxy<CONFIG>>;

template <class CONFIG>
using PyTemplateGalaxy = py::class_<TemplateGalaxy<CONFIG>>;

template <class CONFIG>
using PyTemplateInfo = py::class_<TemplateInfo<CONFIG>>;

template <class FP>
using PyUniformDeviate = py::class_<ran::UniformDeviate<FP>>;

template <class CONFIG>
using PyKDTreePrior = py::class_<KDTreePrior<CONFIG>>;

template <bool FIX_CENTER, bool USE_CONC, bool USE_MAG, int N_COLORS,
          bool USE_FLOAT>
static void declareBfdConfig(py::module &mod, bool fix_center, bool use_conc,
                             bool use_mag, int n_colors, bool use_float) {
  std::string label = optionsToString("BFDConfig", FIX_CENTER, USE_CONC,
                                      USE_MAG, N_COLORS, USE_FLOAT);

  PyBfdConfig<FIX_CENTER, USE_CONC, USE_MAG, N_COLORS, USE_FLOAT> cls(
      mod, label.c_str());
  typedef BfdConfig<FIX_CENTER, USE_CONC, USE_MAG, N_COLORS, USE_FLOAT> BC;
  cls.def(py::init<>());
  cls.attr("MF") = py::int_(BC::MF);
  cls.attr("MR") = py::int_(BC::MR);
  cls.attr("M1") = py::int_(BC::M1);
  cls.attr("M2") = py::int_(BC::M2);
  cls.attr("MC") = py::int_(BC::MC);
  cls.attr("MC0") = py::int_(BC::MC0);
  cls.attr("MSIZE") = py::int_(BC::MSIZE);
  cls.attr("MX") = py::int_(BC::MX);
  cls.attr("MY") = py::int_(BC::MY);
  cls.attr("XYSIZE") = py::int_(BC::XYSIZE);
  cls.attr("MXYSIZE") = py::int_(BC::MXYSIZE);
  cls.attr("G1") = py::int_(BC::G1);
  cls.attr("G2") = py::int_(BC::G2);
  cls.attr("MU") = py::int_(BC::MU);
}

template <class BC>
std::string formatMomentStr(Moment<BC> const &self) {
  std::ostringstream os;
  os << "Even [";
  for (auto i = 0; i < self.m.size(); ++i) os << self.m[i] << ",";
  os.seekp(-1, std::ios_base::end);
  os << "], Odd [";
  for (auto i = 0; i < self.xy.size(); ++i) os << self.xy[i] << ",";
  os.seekp(-1, std::ios_base::end);
  os << "]";
  return os.str();
}

template <bool FIX_CENTER, bool USE_CONC, bool USE_MAG, int N_COLORS,
          bool USE_FLOAT>
static void declareMoment(py::module &mod, bool fix_center, bool use_conc,
                          bool use_mag, int n_colors, bool use_float) {
  std::string label = optionsToString("Moment", FIX_CENTER, USE_CONC, USE_MAG,
                                      N_COLORS, USE_FLOAT);

  typedef BfdConfig<FIX_CENTER, USE_CONC, USE_MAG, N_COLORS, USE_FLOAT> BC;
  typedef typename BC::MVector MVector;
  typedef typename BC::PqrVector PqrVector;
  typedef typename BC::XYVector XYVector;
  typedef typename BC::FP FP;

  PyMoment<BC> cls(mod, label.c_str());
  cls.def(py::init<const MVector &, const XYVector &>(), "m"_a = MVector(FP(0)),
          "xy"_a = XYVector(FP(0)));
  cls.def_readwrite("m", &Moment<BC>::m);
  cls.def_readwrite("xy", &Moment<BC>::xy);
  cls.def("rotate", &Moment<BC>::rotate, "theta"_a);
  cls.def("yflip", &Moment<BC>::yflip);
  cls.def(py::self + py::self);
  cls.def(py::self += py::self);
  cls.def(py::self - py::self);
  cls.def(py::self -= py::self);
  cls.def(py::self * FP());
  cls.def(py::self *= FP());
  cls.def(py::self / FP());
  cls.def(py::self /= FP());
  cls.def("__repr__",
          [](Moment<BC> const &self) { return formatMomentStr(self); });
}

template <class BC>
std::string formatMomentCovStr(MomentCovariance<BC> const &self) {
  std::ostringstream os;
  os << "Even: [";
  for (auto i = 0; i < self.m.nrows(); ++i) {
    if (i > 0) os << "      [";
    for (auto j = 0; j < self.m.cols(); ++j) {
      os << self.m(i, j) << ",";
    }
    os.seekp(-1, std::ios_base::end);
    os << "],";
    if (i == self.m.nrows()) os.seekp(-1, std::ios_base::end);
    os << "\n";
  }

  os << "Odd: [";
  for (auto i = 0; i < self.xy.nrows(); ++i) {
    if (i > 0) os << "     [";
    for (auto j = 0; j < self.xy.cols(); ++j) {
      os << self.xy(i, j) << ",";
    }
    os.seekp(-1, std::ios_base::end);
    os << "],";
    if (i == self.xy.nrows()) os.seekp(-1, std::ios_base::end);
    os << "\n";
  }

  return os.str();
}

template <bool FIX_CENTER, bool USE_CONC, bool USE_MAG, int N_COLORS,
          bool USE_FLOAT>
static void declareMomentCov(py::module &mod, bool fix_center, bool use_conc,
                             bool use_mag, int n_colors, bool use_float) {
  std::string label = optionsToString("MomentCov", FIX_CENTER, USE_CONC,
                                      USE_MAG, N_COLORS, USE_FLOAT);

  typedef BfdConfig<FIX_CENTER, USE_CONC, USE_MAG, N_COLORS, USE_FLOAT> BC;
  typedef typename BC::MMatrix MMatrix;
  typedef typename BC::XYMatrix XYMatrix;
  typedef typename BC::FP FP;

  PyMomentCov<BC> cls(mod, label.c_str());
  cls.def(py::init<>());
  cls.def(py::init<const MMatrix &, const XYMatrix &>(), "m"_a, "xy"_a);
  cls.def_readwrite("m", &MomentCovariance<BC>::m);
  cls.def_readwrite("xy", &MomentCovariance<BC>::xy);
  cls.def("rotate", &MomentCovariance<BC>::rotate, "theta"_a);
  cls.def("yflip", &MomentCovariance<BC>::yflip);
  cls.def(py::self + py::self);
  cls.def(py::self += py::self);
  cls.def(py::self - py::self);
  cls.def(py::self -= py::self);
  cls.def(py::self * FP());
  cls.def(py::self *= FP());
  cls.def(py::self / FP());
  cls.def(py::self /= FP());
  cls.def("__repr__", [](MomentCovariance<BC> const &self) {
    return formatMomentCovStr(self);
  });
}

template <class FP>
static void declareKSigmaWeight(py::module &mod, bool use_float) {
  std::string label = "KSigmaWeight";
  if (use_float)
    label += "F";
  else
    label += "D";

  PyKSigmaWeight<FP> cls(mod, label.c_str());
  cls.def(py::init<FP, int>(), "sigma"_a, "n"_a);
  cls.def("kMax", &KSigmaWeight<FP>::kMax);
  cls.def("__call__", &KSigmaWeight<FP>::operator(), "ksq"_a);
}

template <bool FIX_CENTER, bool USE_CONC, bool USE_MAG, int N_COLORS,
          bool USE_FLOAT>
static void declareKGalaxy(py::module &mod, bool fix_center, bool use_conc,
                           bool use_mag, int n_colors, bool use_float) {
  std::string label = optionsToString("KGalaxy", FIX_CENTER, USE_CONC, USE_MAG,
                                      N_COLORS, USE_FLOAT);

  typedef BfdConfig<FIX_CENTER, USE_CONC, USE_MAG, N_COLORS, USE_FLOAT> BC;
  typedef typename BC::FP FP;
  typedef linalg::Vector<FP> Vector;
  typedef linalg::Vector<std::complex<FP>> CVector;
  typedef linalg::Matrix<FP> Matrix;
  typedef linalg::Matrix<std::complex<FP>> CMatrix;
  typedef linalg::DVector2 DVector2;

  PyKGalaxy<BC> cls(mod, label.c_str());
  cls.def(py::init<const KSigmaWeight<FP> &, const CVector &, const Vector &,
                   const Vector &, FP, const std::set<int>, const DVector2 &,
                   long>(),
          "kw"_a, "kval"_a, "kx"_a, "ky"_a, "d2k"_a, "unconjugated"_a,
          "posn_"_a = DVector2(0.), "id_"_a = 0L);
  cls.def(py::init<const KSigmaWeight<FP> &, const CVector &, const Vector &,
                   const Vector &, const Vector &, FP, const std::set<int>,
                   const DVector2 &, const long>(),
          "kw"_a, "kval"_a, "kx"_a, "ky"_a, "kvar_"_a, "d2k"_a,
          "unconjugated"_a, "posn_"_a = DVector2(0.), "id_"_a = 0L);
  cls.def("getSheared", &KGalaxy<BC>::getSheared, "g1"_a, "g2"_a, "mu"_a);
  cls.def("getMoment", &KGalaxy<BC>::getMoment);
  cls.def("getShifted", &KGalaxy<BC>::getShifted, "dx"_a, "dy"_a);
  cls.def("getNulled", [](const KGalaxy<BC> &self, FP maxShift) {
    double dx, dy;
    KGalaxy<BC> *gal = self.getNulled(dx, dy, maxShift);
    return std::make_tuple(gal, dx, dy);
  });
  cls.def("getNulled", &KGalaxy<BC>::getNulled, "dx"_a, "dy"_a, "maxShift"_a);
  cls.def("getTarget", &KGalaxy<BC>::getTarget, "fillCovariance"_a = true);
  cls.def("getTemplate", &KGalaxy<BC>::getTemplate);
  cls.def("getShearedTarget", &KGalaxy<BC>::getShearedTarget, "g1"_a, "g2"_a,
          "mu"_a, "fillCovariance"_a = true);
}

template <bool FIX_CENTER, bool USE_CONC, bool USE_MAG, int N_COLORS,
          bool USE_FLOAT>
static void declareTargetGalaxy(py::module &mod, bool fix_center, bool use_conc,
                                bool use_mag, int n_colors, bool use_float) {
  std::string label = optionsToString("TargetGalaxy", FIX_CENTER, USE_CONC,
                                      USE_MAG, N_COLORS, USE_FLOAT);

  typedef BfdConfig<FIX_CENTER, USE_CONC, USE_MAG, N_COLORS, USE_FLOAT> BC;
  typedef typename BC::FP FP;
  typedef linalg::Vector<FP> Vector;
  typedef linalg::Vector<std::complex<FP>> CVector;
  typedef linalg::Matrix<FP> Matrix;
  typedef linalg::Matrix<std::complex<FP>> CMatrix;
  typedef linalg::DVector DVector;
  typedef linalg::DVector2 DVector2;
  typedef TemplateGalaxy<BC> TG;

  PyTargetGalaxy<BC> cls(mod, label.c_str());
  cls.def(py::init<const Moment<BC> &, const MomentCovariance<BC> &,
                   const DVector &, const long>(),
          "mom_"_a = Moment<BC>(), "cov_"_a = MomentCovariance<BC>(),
          "position"_a = DVector2(0.), "id_"_a = 0L);
  cls.def_readwrite("mom", &TargetGalaxy<BC>::mom);
  cls.def_readwrite("cov", &TargetGalaxy<BC>::cov);
  cls.def_readwrite("id", &TargetGalaxy<BC>::id);
  cls.def_readwrite("position", &TargetGalaxy<BC>::position);
  cls.def("rotate", &TargetGalaxy<BC>::rotate, "theta"_a);
  cls.def("yflip", &TargetGalaxy<BC>::yflip);
}

template <bool FIX_CENTER, bool USE_CONC, bool USE_MAG, int N_COLORS,
          bool USE_FLOAT>
static void declareTemplateGalaxy(py::module &mod, bool fix_center,
                                  bool use_conc, bool use_mag, int n_colors,
                                  bool use_float) {
  std::string label = optionsToString("TemplateGalaxy", FIX_CENTER, USE_CONC,
                                      USE_MAG, N_COLORS, USE_FLOAT);

  typedef BfdConfig<FIX_CENTER, USE_CONC, USE_MAG, N_COLORS, USE_FLOAT> BC;
  typedef typename BC::FP FP;
  typedef std::complex<FP> FPC;
  typedef linalg::Vector<FP> Vector;
  typedef linalg::Vector<std::complex<FP>> CVector;
  typedef linalg::Matrix<FP> Matrix;
  typedef linalg::Matrix<std::complex<FP>> CMatrix;
  typedef linalg::DVector DVector;
  typedef linalg::Matrix<std::complex<FP>> DerivMatrix;
  typedef TemplateGalaxy<BC> TG;

  PyTemplateGalaxy<BC> cls(mod, label.c_str());
  cls.def(py::init<const DerivMatrix &, const DerivMatrix &, const FP,
                   const long, const FP>(),
          "mderiv_"_a = DerivMatrix(BC::MSIZE, BC::DSIZE, 0.),
          "xyderiv_"_a = DerivMatrix(BC::XYSIZE, BC::DSIZE, 0.), "nda_"_a = 1.,
          "id_"_a = 0L, "jSupression"_a = 1.);
  cls.def_readwrite("mDeriv", &TemplateGalaxy<BC>::mDeriv);
  cls.def_readwrite("xyDeriv", &TemplateGalaxy<BC>::xyDeriv);
  cls.def_readwrite("nda", &TemplateGalaxy<BC>::nda);
  cls.def_readwrite("id", &TemplateGalaxy<BC>::id);
  cls.def_readwrite("jSuppression", &TemplateGalaxy<BC>::jSuppression);
  cls.def("rotate", &TemplateGalaxy<BC>::rotate, "theta"_a);
  cls.def("yflip", &TemplateGalaxy<BC>::yflip);
  cls.def("realMDerivs", &TemplateGalaxy<BC>::realMDerivs);
  cls.def("realXYDerivs", &TemplateGalaxy<BC>::realXYDerivs);
  cls.attr("MF") = py::int_(TG::MF);
  cls.attr("MR") = py::int_(TG::MR);
  cls.attr("ME") = py::int_(TG::ME);
  cls.attr("MC") = py::int_(TG::MC);
  cls.attr("MC0") = py::int_(TG::MC0);
  cls.attr("MSIZE") = py::int_(TG::MSIZE);
  cls.attr("MX") = py::int_(TG::MX);
  cls.attr("XYSIZE") = py::int_(TG::XYSIZE);
  cls.attr("D0") = py::int_(TG::D0);
  cls.attr("Q") = py::int_(TG::Q);
  cls.attr("DU") = py::int_(TG::DU);
  cls.attr("DVb") = py::int_(TG::DVb);
  cls.attr("R") = py::int_(TG::R);
  cls.attr("DU_DU") = py::int_(TG::DU_DU);
  cls.attr("DU_DV") = py::int_(TG::DU_DV);
  cls.attr("DU_DVb") = py::int_(TG::DU_DVb);
  cls.attr("DV_DV") = py::int_(TG::DV_DV);
  cls.attr("DVb_DVb") = py::int_(TG::DVb_DVb);
  cls.attr("DV_DVb") = py::int_(TG::DVb_DVb);
  cls.attr("DSIZE") = py::int_(TG::DSIZE);
}


template <bool FIX_CENTER, bool USE_CONC, bool USE_MAG, int N_COLORS,
          bool USE_FLOAT>
static void declareTemplateInfo(py::module &mod, bool fix_center,
                                bool use_conc, bool use_mag, int n_colors,
                                bool use_float) {
  std::string label = optionsToString("TemplateInfo", FIX_CENTER, USE_CONC,
                                      USE_MAG, N_COLORS, USE_FLOAT);

  typedef BfdConfig<FIX_CENTER, USE_CONC, USE_MAG, N_COLORS, USE_FLOAT> BC;

  PyTemplateInfo<BC> cls(mod, label.c_str());
  cls.def(py::init<>());
  cls.def(py::init<const TemplateGalaxy<BC>& >(), "tmpl"_a);
  cls.def_readwrite("dm", &TemplateInfo<BC>::dm);
  cls.def_readwrite("dxy", &TemplateInfo<BC>::dxy);
  cls.def_readwrite("nda", &TemplateInfo<BC>::nda);
  cls.def_readwrite("id", &TemplateInfo<BC>::id);
  cls.def_readwrite("m", &TemplateInfo<BC>::m);
  cls.def("getM", &TemplateInfo<BC>::getM);
  cls.def("setM", &TemplateInfo<BC>::setM, "_m"_a);
}

template <class FP>
static void declareUniformDeviate(py::module &mod, std::string label) {
  PyUniformDeviate<FP> cls(mod, label.c_str());
  cls.def(py::init<>());
  cls.def(py::init<const long>(), "lseed"_a);
  cls.def("__call__", &ran::UniformDeviate<FP>::operator());
}

template <class CONFIG>
class PqrWrapper {
 public:
  PqrWrapper() {}
  PqrWrapper(Pqr<CONFIG> _p) : p(_p) {}
  Pqr<CONFIG> p;
};

template <bool FIX_CENTER, bool USE_CONC, bool USE_MAG, int N_COLORS,
          bool USE_FLOAT>
static void declareKDTreePrior(py::module &mod, bool fix_center, bool use_conc,
                               bool use_mag, int n_colors, bool use_float) {
  std::string label = optionsToString("KDTreePrior", FIX_CENTER, USE_CONC,
                                      USE_MAG, N_COLORS, USE_FLOAT);

  typedef BfdConfig<FIX_CENTER, USE_CONC, USE_MAG, N_COLORS, USE_FLOAT> BC;
  typedef typename BC::FP FP;
  typedef std::complex<FP> FPC;
  typedef linalg::Vector<FP> Vector;
  typedef linalg::Vector<std::complex<FP>> CVector;
  typedef linalg::Matrix<FP> Matrix;
  typedef linalg::Matrix<std::complex<FP>> CMatrix;
  typedef linalg::DVector DVector;
  typedef linalg::Matrix<std::complex<FP>> DerivMatrix;

  PyKDTreePrior<BC> cls(mod, label.c_str());
  cls.def(py::init<FP, FP, const MomentCovariance<BC> &,
                   ran::UniformDeviate<double> &, int, bool, FP, FP, FP, FP,
                   bool, bool>(),
          "fluxMin_"_a, "fluxMax_"_a, "nominalCov"_a, "ud_"_a, "nSamp_"_a,
          "selectionOnly_"_a = false, "noiseFactor"_a = 1., "sigmaStep_"_a = 1.,
          "sigmaCutoff_"_a = 6.5, "sigmaBuffer_"_a = 1.,
          "invariantCovariance_"_a = false, "fixedNoise_"_a = false);
  cls.def("prepare", &KDTreePrior<BC>::prepare);
  // cls.def("getPqr", &KDTreePrior<BC>::getPqr, "gal"_a, "nTemplates"_a,
  // "nUnique"_a) ;
  cls.def("getPqr", [](const KDTreePrior<BC> &self, const TargetGalaxy<BC> &gal,
                       int &nTemplates, int &nUnique) {
    Pqr<BC> pqr = self.getPqr(gal, nTemplates, nUnique);
    return std::make_tuple(PqrWrapper<BC>(pqr), nTemplates, nUnique);
  });
  cls.def("setSampleWeights", &KDTreePrior<BC>::setSampleWeights, "b"_a);
}

template <class BC>
std::string formatPqrStr(PqrWrapper<BC> const &self) {
  std::ostringstream os;
  os << "[";
  for (auto i = 0; i < self.p.size(); ++i) os << self.p[i] << ",";
  os.seekp(-1, std::ios_base::end);
  os << "]";
  return os.str();
}

template <bool FIX_CENTER, bool USE_CONC, bool USE_MAG, int N_COLORS,
          bool USE_FLOAT>
static void declarePqr(py::module &mod, bool fix_center, bool use_conc,
                       bool use_mag, int n_colors, bool use_float) {
  std::string label = optionsToString("Pqr", FIX_CENTER, USE_CONC, USE_MAG,
                                      N_COLORS, USE_FLOAT);

  typedef BfdConfig<FIX_CENTER, USE_CONC, USE_MAG, N_COLORS, USE_FLOAT> BC;
  typedef typename BC::FP FP;
  typedef typename BC::PqrVector PqrVector;
  typedef typename BC::QVector QVector;
  typedef typename BC::RMatrix RMatrix;
  typedef linalg::DVector DVector;

  py::class_<PqrWrapper<BC>> cls(mod, label.c_str());
  cls.def(py::init<PqrVector>(), "pqr"_a = PqrVector());
  cls.def_readwrite("_pqr", &PqrWrapper<BC>::p);
  cls.def("getP", [](PqrWrapper<BC> const &self) { return self.p.getP(); });
  cls.def("getQ", [](PqrWrapper<BC> const &self) { return self.p.getQ(); });
  cls.def("getR", [](PqrWrapper<BC> const &self) { return self.p.getR(); });
  cls.def("setP", [](PqrWrapper<BC> &self, FP &val) { self.p.setP(val); });
  cls.def("setQ", [](PqrWrapper<BC> &self, QVector &val) { self.p.setQ(val); });
  cls.def("chainRule", [](PqrWrapper<BC> &self, FP f, FP df, FP ddf) {
    return PqrWrapper<BC>(self.p.chainRule(f, df, ddf));
  });
  cls.def("__mul__",
          [](PqrWrapper<BC> &self, const PqrWrapper<BC> &other) {
            return PqrWrapper<BC>(self.p *= other.p);
          },
          py::is_operator());
  cls.def("neglog",
          [](PqrWrapper<BC> &self) { PqrWrapper<BC>(self.p.neglog()); });
  cls.def("rotate",
          [](PqrWrapper<BC> &self, double theta) { self.p.rotate(theta); });
  cls.def("getG", [](PqrWrapper<BC> &self) {
    QVector g;
    RMatrix cov;
    self.p.getG(g, cov);
    return std::make_tuple(g, cov);
  });
  cls.def("__call__",
          [](PqrWrapper<BC> &self, QVector &q) { return self.p(q); });
  cls.def("__repr__",
          [](PqrWrapper<BC> const &self) { return formatPqrStr(self); });
}

PYBIND11_MODULE(pybind_bfd, mod) {
  declareKSigmaWeight<float>(mod, true);
  declareUniformDeviate<double>(mod, "UniformDeviate");

  // Fix-center=False, concentration=False, mag=False, colors=0, float=True
  declareBfdConfig<false, false, false, 0, true>(mod, false, false, false, 0,
                                                 true);
  declareMoment<false, false, false, 0, true>(mod, false, false, false, 0,
                                              true);
  declareMomentCov<false, false, false, 0, true>(mod, false, false, false, 0,
                                                 true);
  declareKGalaxy<false, false, false, 0, true>(mod, false, false, false, 0,
                                               true);
  declareTargetGalaxy<false, false, false, 0, true>(mod, false, false, false, 0,
                                                    true);
  declareTemplateGalaxy<false, false, false, 0, true>(mod, false, false, false,
                                                      0, true);
  declareTemplateInfo<false, false, false, 0, true>(mod, false, false, false,
                                                      0, true);
  declareKDTreePrior<false, false, false, 0, true>(mod, false, false, false, 0,
                                                   true);
  declarePqr<false, false, false, 0, true>(mod, false, false, false, 0, true);


  // Fix-center=False, concentration=False, mag=True, colors=0, float=True
  declareBfdConfig<false, false, true, 0, true>(mod, false, false, true, 0,
                                                true);
  declareMoment<false, false, true, 0, true>(mod, false, false, true, 0, true);
  declareMomentCov<false, false, true, 0, true>(mod, false, false, true, 0,
                                                true);
  declareKGalaxy<false, false, true, 0, true>(mod, false, false, true, 0, true);
  declareTargetGalaxy<false, false, true, 0, true>(mod, false, false, true, 0,
                                                   true);
  declareTemplateGalaxy<false, false, true, 0, true>(mod, false, false, true, 0,
                                                     true);
  declareTemplateInfo<false, false, true, 0, true>(mod, false, false, true, 0,
                                                     true);
  declareKDTreePrior<false, false, true, 0, true>(mod, false, false, true, 0,
                                                  true);
  declarePqr<false, false, true, 0, true>(mod, false, false, true, 0, true);

  // Fix-center=False, concentration=True, mag=True, colors=0, float=True
  declareBfdConfig<false, true, true, 0, true>(mod, false, true, true, 0, true);
  declareMoment<false, true, true, 0, true>(mod, false, true, true, 0, true);
  declareMomentCov<false, true, true, 0, true>(mod, false, true, true, 0, true);
  declareKGalaxy<false, true, true, 0, true>(mod, false, true, true, 0, true);
  declareTargetGalaxy<false, true, true, 0, true>(mod, false, true, true, 0,
                                                  true);
  declareTemplateGalaxy<false, true, true, 0, true>(mod, false, true, true, 0,
                                                    true);
  declareTemplateInfo<false, true, true, 0, true>(mod, false, true, true, 0,
                                                    true);
  declareKDTreePrior<false, true, true, 0, true>(mod, false, true, true, 0,
                                                 true);
  declarePqr<false, true, true, 0, true>(mod, false, true, true, 0, true);

  // Add 5 color fluxes for use in LSST
  declareBfdConfig<false, true, true, 5, true>(mod, false, true, true, 5, true);
  declareMoment<false, true, true, 5, true>(mod, false, true, true, 5, true);
  declareMomentCov<false, true, true, 5, true>(mod, false, true, true, 5, true);
  declareKGalaxy<false, true, true, 5, true>(mod, false, true, true, 5, true);
  declareTargetGalaxy<false, true, true, 5, true>(mod, false, true, true, 5,
                                                  true);
  declareTemplateGalaxy<false, true, true, 5, true>(mod, false, true, true, 5,
                                                    true);
  declareTemplateInfo<false, true, true, 5, true>(mod, false, true, true, 5,
                                                    true);
  declareKDTreePrior<false, true, true, 5, true>(mod, false, true, true, 5,
                                                 true);
  declarePqr<false, true, true, 5, true>(mod, false, true, true, 5, true);
}

}  // namespace

}  // namespace bfd
}  // namespace desc
}  // namespace lsst