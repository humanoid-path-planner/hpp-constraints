// Copyright (c) 2017, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

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

        bool isEqual(const DifferentiableFunction& other) const {
          const ActiveSetDifferentiableFunction& castother = dynamic_cast<const ActiveSetDifferentiableFunction&>(other);
          if (!DifferentiableFunction::isEqual(other))
            return false;
          
          if (function_ != castother.function_)
            return false;
          if (intervals_ != castother.intervals_)
            return false;
          
          return true;
        }

        DifferentiableFunctionPtr_t function_;
        segments_t intervals_;
    }; // class ActiveSetDifferentiableFunction
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_ACTIVE_SET_DIFFERENTIABLE_FUNCTION_HH
