// Copyright (c) 2014, LAAS-CNRS
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//
// This file is part of hpp-constraints.
// hpp-constraints is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-constraints is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-constraints. If not, see <http://www.gnu.org/licenses/>.

#ifndef HPP_CONSTRAINTS_TOOL_HH
#define HPP_CONSTRAINTS_TOOL_HH

#include "hpp/constraints/fwd.hh"

namespace hpp {
  namespace constraints {
    void computeLog (vectorOut_t result, double& theta, const fcl::Matrix3f& Rerror);

    void computeJlog (const double& theta, vectorIn_t r, eigen::matrix3_t& Jlog);

    template < typename VectorType, typename MatrixType >
    static void computeCrossMatrix (const VectorType& v, MatrixType& m)
    {
      m.diagonal ().setZero ();
      m (0,1) = -v [2]; m (1,0) = v [2];
      m (0,2) = v [1]; m (2,0) = -v [1];
      m (1,2) = -v [0]; m (2,1) = v [0];
    }
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_TOOL_HH
