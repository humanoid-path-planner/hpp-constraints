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

#ifndef HPP_CONSTRAINTS_FUNCTION_DIFFERENCE_HH
# define HPP_CONSTRAINTS_FUNCTION_DIFFERENCE_HH

# include <hpp/constraints/config.hh>
# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/fwd.hh>

namespace hpp {
  namespace constraints {
    namespace function {
    class Difference;
    typedef boost::shared_ptr<Difference> DifferencePtr_t;

    /// Compute the difference between the value of the function in two points.
    /// i.e.: \f$ f (q_0, ... , q_n) = f_{inner} (q_{left}) - f_{inner} (q_{right}) \f$
    class HPP_CONSTRAINTS_LOCAL Difference :
      public constraints::DifferentiableFunction
    {
      public:
        typedef boost::shared_ptr<Difference> Ptr_t;
        Difference (const DifferentiableFunctionPtr_t& inner,
            const size_type& nArgs, const size_type& nDers,
            const segment_t& lInArgs, const segment_t& lInDers,
            const segment_t& rInArgs, const segment_t& rInDers);

      protected:
        void impl_compute (LiegroupElementRef y, vectorIn_t arg) const
        {
          inner_->value(l_, arg.segment (lsa_.first, lsa_.second));
          inner_->value(r_, arg.segment (rsa_.first, rsa_.second));
          y.vector() = l_ - r_;
        }

        void impl_jacobian (matrixOut_t J, vectorIn_t arg) const
        {
          inner_->jacobian(
              J.middleCols (lsd_.first, lsd_.second),
              arg.segment (lsa_.first, lsa_.second));
          inner_->jacobian(
              J.middleCols (rsd_.first, rsd_.second),
              arg.segment (rsa_.first, rsa_.second));
          J.middleCols (rsd_.first, rsd_.second) *= -1;
        }

        std::ostream& print (std::ostream& os) const;

        DifferentiableFunctionPtr_t inner_;
        const segment_t lsa_, lsd_;
        const segment_t rsa_, rsd_;

        mutable LiegroupElement l_, r_;
    }; // class Difference
    } // namespace function
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_FUNCTION_DIFFERENCE_HH
