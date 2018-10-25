// Copyright (c) 2018, Joseph Mirabel
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

#include <hpp/constraints/matrix-view.hh>
#include <hpp/constraints/function/difference.hh>

namespace hpp {
  namespace constraints {
    namespace function {
      namespace {
        std::string toStr (const segment_t& s)
        {
          std::ostringstream os;
          os << pretty_print (s);
          return os.str();
        }
      }

      Difference::Difference (const DifferentiableFunctionPtr_t& inner,
          const size_type& nArgs, const size_type& nDers,
          const segment_t& lInArgs, const segment_t& lInDers,
          const segment_t& rInArgs, const segment_t& rInDers) :
        DifferentiableFunction (nArgs, nDers, LiegroupSpace::Rn(inner->outputSpace()->nv()), inner->name() + " | " + toStr(lInArgs) + " - " + toStr(rInArgs)),
        inner_ (inner),
        lsa_ (lInArgs), lsd_ (lInDers),
        rsa_ (rInArgs), rsd_ (rInDers),
        l_ (inner->outputSpace()), r_ (inner->outputSpace())
      {
        activeParameters_.setConstant(false);
        activeParameters_.segment(lsa_.first, lsa_.second)
          = inner_->activeParameters();
        activeParameters_.segment(rsa_.first, rsa_.second)
          = inner_->activeParameters();

        activeDerivativeParameters_.setConstant(false);
        activeDerivativeParameters_.segment(lsd_.first, lsd_.second)
          = inner_->activeDerivativeParameters();
        activeDerivativeParameters_.segment(rsd_.first, rsd_.second)
          = inner_->activeDerivativeParameters();
      }

      std::ostream& Difference::print (std::ostream& os) const
      {
        constraints::DifferentiableFunction::print(os);
        return os << incindent << iendl << *inner_ << decindent;
      }
    } // namespace function
  } // namespace constraints
} // namespace hpp
