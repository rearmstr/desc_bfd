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

#ifndef LSST_MEAS_EXTENSION_BFD_bfd_h_INCLUDED
#define LSST_MEAS_EXTENSION_BFD_bfd_h_INCLUDED

#include "Eigen/Core"
#include "ndarray_fwd.h"
#include "lsst/afw/table/fwd.h"

namespace lsst { namespace desc { namespace bfd {

// This corresponds to all 6 moments
const int UseMoments = 15;
typedef float Pixel;
typedef float Scalar;
const int NMoment = 6;
const int Ng = 2;
const int Npqr = 6;
typedef Eigen::Matrix<Scalar,NMoment,NMoment> Matrix;
typedef Eigen::Matrix<Scalar,NMoment,1> Vector;
typedef afw::table::Key<Scalar> ScalarKey;
typedef afw::table::Key< afw::table::Array<Scalar> > ArrayKey;

}}} // namespace lsst::desc::bfd

#endif // !LSST_MEAS_EXTENSION_BFD_bfd_h_INCLUDED
