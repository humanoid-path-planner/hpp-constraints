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

#ifndef HPP_CONSTRAINTS_SOLVER_HH
#define HPP_CONSTRAINTS_SOLVER_HH

#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/config.hh>

namespace hpp {
  namespace constraints {
    class HPP_CONSTRAINTS_DLLAPI Solver
    {
      public:
        enum Type {
          Explicit,
          Iterative
        };

        virtual bool solve (vectorOut_t arg) const = 0;

        const Type& type () const
        {
          return type_;
        }

      protected:
        Solver (const Type& t) : type_ (t) {}

        const Type type_;
    };
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_SOLVER_HH
