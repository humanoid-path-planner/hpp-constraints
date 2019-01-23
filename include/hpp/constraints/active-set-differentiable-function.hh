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

#ifndef HPP_CONSTRAINTS_ACTIVE_SET_DIFFERENTIABLE_FUNCTION_HH
#define HPP_CONSTRAINTS_ACTIVE_SET_DIFFERENTIABLE_FUNCTION_HH

#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/config.hh>
#include <hpp/constraints/differentiable-function.hh>

namespace hpp {
  namespace constraints {

    /// Handle bounds on input variables of a differentiable function.
    ///
    /// This class is a decorator of class DifferentiableFunction that
    /// sets to 0 some columns of the Jacobian of the function.
    ///
    /// The class is used to handle saturation of input variables of
    /// the function during numerical resolution of implicit constraints
    /// built with the function.
    class HPP_CONSTRAINTS_DLLAPI ActiveSetDifferentiableFunction :
      public DifferentiableFunction
    {
      public:
        /// Constructor
        /// \param f initial differentiable function,
        /// \param intervals set of intervals of indices corresponding to saturated
        ///            input variables.
        ActiveSetDifferentiableFunction (const DifferentiableFunctionPtr_t& f,
            segments_t intervals)
          : DifferentiableFunction(
              f->inputSize(), f->inputDerivativeSize(),
              f->outputSpace (),
              "ActiveSet_on_" + f->name ())
          , function_(f)
          , intervals_(intervals)
        {
          context (f->context());
        }

        /// Get the original function
        const DifferentiableFunction& function() const
        {
          return *function_;
        }

        /// Get the original function
        const DifferentiableFunctionPtr_t& functionPtr() const
        {
          return function_;
        }

      protected:
        typedef std::vector < segments_t > intervalss_t;

        /// User implementation of function evaluation
        virtual void impl_compute (LiegroupElementRef result,
                                   vectorIn_t argument) const
        {
          function_->value(result, argument);
        }

        virtual void impl_jacobian (matrixOut_t jacobian,
                                    vectorIn_t arg) const
        {
          function_->jacobian(jacobian, arg);
          for (segments_t::const_iterator _int = intervals_.begin ();
              _int != intervals_.end (); ++_int)
            jacobian.middleCols (_int->first, _int->second).setZero ();
        }

        DifferentiableFunctionPtr_t function_;
        segments_t intervals_;
    }; // class ActiveSetDifferentiableFunction
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_ACTIVE_SET_DIFFERENTIABLE_FUNCTION_HH
