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
    void computeLog (vectorOut_t result, double& theta, const fcl::Matrix3f& Rerror)
    {
      double tr = Rerror (0, 0) + Rerror (1, 1) + Rerror (2, 2);
      if (tr > 3) tr = 3;
      if (tr < -1) tr = -1;
      theta = acos ((tr - 1)/2);
      if (theta > 1e-6) {
        result [0] = theta*(Rerror (2, 1) - Rerror (1, 2))/(2*sin(theta));
        result [1] = theta*(Rerror (0, 2) - Rerror (2, 0))/(2*sin(theta));
        result [2] = theta*(Rerror (1, 0) - Rerror (0, 1))/(2*sin(theta));
      } else {
        result [0] = (Rerror (2, 1) - Rerror (1, 2))/2;
        result [1] = (Rerror (0, 2) - Rerror (2, 0))/2;
        result [2] = (Rerror (1, 0) - Rerror (0, 1))/2;
      }
    }

    void computeJlog (const double& theta, vectorIn_t r, eigen::matrix3_t& Jlog)
    {
      if (theta < 1e-6)
	Jlog.setIdentity ();
      else {
	Jlog.setZero ();
	// Jlog = alpha I
	double alpha = .5*theta*sin(theta)/(1-cos (theta));
	Jlog (0,0) = Jlog (1,1) = Jlog (2,2) = alpha;
	// Jlog += -r_{\times}/2
	Jlog (0,1) = .5*r [2]; Jlog (1,0) = -.5*r [2];
	Jlog (0,2) = -.5*r [1]; Jlog (2,0) = .5*r [1];
	Jlog (1,2) = .5*r [0]; Jlog (2,1) = -.5*r [0];
	alpha = 1/(theta*theta) - sin(theta)/(2*theta*(1-cos(theta)));
	Jlog += alpha * r * r.transpose ();
      }
    }
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_TOOL_HH
