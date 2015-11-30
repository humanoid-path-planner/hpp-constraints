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

#include <hpp/constraints/differentiable-function.hh>

namespace hpp {
  namespace constraints {
      void DifferentiableFunction::finiteDifferenceForward
        (matrixOut_t jacobian, vectorIn_t x, value_type eps) const
      {
        using std::abs;

        value_type h;
        size_type n = x.size();
        vector_t x_dx = x;
        vector_t f_x   (outputSize()),
                 f_x_dx(outputSize());
        impl_compute (f_x, x);

        for (size_type j = 0; j < n; ++j) {
          h = eps * abs(x[j]);
          if (h == 0) h = eps;
          x_dx[j] += h;
          impl_compute (f_x_dx, x_dx);
          jacobian.col (j) = (f_x_dx - f_x) / h;
          x_dx[j] = x[j];
        }
      }

      void DifferentiableFunction::finiteDifferenceCentral
        (matrixOut_t jacobian, vectorIn_t x, value_type eps) const
      {
        using std::abs;

        value_type h;
        size_type n = x.size();
        vector_t x_dx = x;
        vector_t f_x_mdx (outputSize()),
                 f_x_pdx (outputSize());

        for (size_type j = 0; j < n; ++j) {
          h = eps * abs(x[j]);
          if (h == 0) h = eps;
          x_dx[j] -= h;
          impl_compute (f_x_mdx, x_dx);
          x_dx[j] = x[j] + h;
          impl_compute (f_x_pdx, x_dx);
          jacobian.col (j) = (f_x_pdx - f_x_mdx) / (2*h);
          x_dx[j] = x[j];
        }
      }
  } // namespace constraints
} // namespace hpp
