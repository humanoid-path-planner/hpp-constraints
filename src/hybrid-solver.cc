// Copyright (c) 2017, Joseph Mirabel
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

#include <hpp/constraints/hybrid-solver.hh>

namespace hpp {
  namespace constraints {
    void HybridSolver::explicitSolverHasChanged()
    {
      reduction(explicit_.inDers());
    }

    void HybridSolver::updateJacobian (vectorIn_t arg) const
    {
      // Compute Je_
      explicit_.jacobian(JeExpanded_, arg);
      explicit_.viewJacobian(JeExpanded_).writeTo(Je_);

      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        Data& d = datas_[i];
        d.reducedJ.noalias() += explicit_.outDers().rviewTranspose(d.jacobian).eval() * Je_;
      }
    }
  } // namespace constraints
} // namespace hpp
