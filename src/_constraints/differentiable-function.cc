// Copyright (c) 2015, Joseph Mirabel
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

#include <hpp/_constraints/differentiable-function.hh>

#include <hpp/model/configuration.hh>

namespace hpp {
  namespace _constraints {
      void DifferentiableFunction::finiteDifferenceForward
        (matrixOut_t jacobian, vectorIn_t x,
         DevicePtr_t robot, value_type eps) const
      {
        // TODO: When robot is not null, it would be better to iterate over the
        // joints.
        using std::abs;

        size_type n = inputDerivativeSize();
        vector_t h = vector_t::Zero (inputDerivativeSize ());
        vector_t x_dx = x;
        vector_t f_x   (outputSize()),
                 f_x_dx(outputSize());
        impl_compute (f_x, x);

        for (size_type j = 0; j < n; ++j) {
          if (robot) {
            JointPtr_t jt = robot->getJointAtVelocityRank (j);
            h[j] = eps * x.segment (jt->rankInConfiguration (),
                                    jt->configSize ()).norm();
          }
          else h[j] = eps * abs(x[j]);
          if (h[j] == 0) h[j] = eps;

          if (robot) integrate (robot, x, h, x_dx);
          else x_dx[j] += h[j];

          impl_compute (f_x_dx, x_dx);
          jacobian.col (j) = (f_x_dx - f_x) / h[j];
          if (jacobian.col(j).hasNaN ()) {
            hppDout (error, "Forward finite difference: NaN");
          }
          x_dx[j] = x[j];
          h[j] = 0;
        }
        if (jacobian.hasNaN ()) {
          hppDout (error, "Forward finite difference: NaN");
        }
      }

      void DifferentiableFunction::finiteDifferenceCentral
        (matrixOut_t jacobian, vectorIn_t x,
         DevicePtr_t robot, value_type eps) const
      {
        using std::abs;
        using hpp::model::integrate;

        size_type n = inputDerivativeSize();
        vector_t x_dx = x;
        vector_t h = vector_t::Zero (inputDerivativeSize ());
        vector_t f_x_mdx (outputSize()),
                 f_x_pdx (outputSize());

        for (size_type j = 0; j < n; ++j) {
          if (robot) {
            JointPtr_t jt = robot->getJointAtVelocityRank (j);
            h[j] = eps * x.segment (jt->rankInConfiguration (),
                                    jt->configSize ()).norm();
          }
          else h[j] = eps * abs(x[j]);
          if (h[j] == 0) h[j] = eps;

          if (robot) integrate (robot, x, -h, x_dx);
          else x_dx[j] -= h[j];
          impl_compute (f_x_mdx, x_dx);

          if (robot) integrate (robot, x, h, x_dx);
          else x_dx[j] = x[j] + h[j];
          impl_compute (f_x_pdx, x_dx);

          jacobian.col (j) = (f_x_pdx - f_x_mdx) / (2*h[j]);
          x_dx[j] = x[j];
          h[j] = 0;
        }
        if (jacobian.hasNaN ()) {
          hppDout (error, "Central finite difference: NaN");
        }
      }
  } // namespace _constraints
} // namespace hpp
