// Copyright (c) 2017 - 2018, CNRS
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr), Florent Lamiraux
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

#include <hpp/constraints/function/of-parameter-subset.hh>

#include <hpp/util/indent.hh>

namespace hpp {
  namespace constraints {
    namespace function {
      namespace {
        std::string toStr (const segment_t& s)
        {
          std::ostringstream os;
          os << "[ " << s.first << ", " << s.first + s.second << " ]";
          return os.str();
        }
      }

      OfParameterSubset::OfParameterSubset (const DifferentiableFunctionPtr_t& g,
                          const size_type& nArgs, const size_type& nDers,
                          const segment_t& inArgs, const segment_t& inDers) :
        DifferentiableFunction (nArgs, nDers, g->outputSpace(),
                                g->name() + " | " + toStr(inArgs)),
        g_ (g), sa_ (inArgs), sd_ (inDers)
      {
        activeParameters_.setConstant(false);
        activeParameters_.segment(sa_.first, sa_.second)
          = g_->activeParameters();

        activeDerivativeParameters_.setConstant(false);
        activeDerivativeParameters_.segment(sd_.first, sd_.second)
          = g_->activeDerivativeParameters();
      }

      std::ostream& OfParameterSubset::print (std::ostream& os) const
      {
        constraints::DifferentiableFunction::print(os);
        return os << incindent << iendl << *g_ << decindent;
      }
    } // namespace function
  } // namespace constraints
} // namespace hpp
