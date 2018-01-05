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

#include <hpp/constraints/differentiable-function-stack.hh>

#include <hpp/util/indent.hh>

namespace hpp {
  namespace constraints {
    std::ostream& DifferentiableFunctionStack::print (std::ostream& os) const
    {
      DifferentiableFunction::print (os) << incindent;
      for (std::size_t i = 0; i < functions_.size(); ++i)
        os << iendl << i << ": " << *functions_[i];
      return os << decindent;
    }
  } // namespace constraints
} // namespace hpp
