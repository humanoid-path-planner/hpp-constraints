// Copyright (c) 2015, LAAS-CNRS
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

#include "hpp/constraints/implicit.hh"
#include <hpp/constraints/differentiable-function.hh>

namespace hpp {
  namespace constraints {
    size_type computeParameterSize (const ComparisonTypes_t& comp)
    {
      size_type size = 0;
      for (std::size_t i = 0; i < comp.size(); ++i) {
        switch (comp[i]) {
        case Superior:
        case Inferior:
          break;
        case Equality:
          ++size;
          break;
        default:
          break;
        }
      }
      return size;
    }
    void Implicit::rightHandSide (vectorIn_t rightHandSide)
    {
      LiegroupSpacePtr_t space (function_->outputSpace());
      assert (rightHandSide.size () == space->nq ());
      pinocchio::LiegroupElementConstRef output
        (space->elementConstRef (rightHandSide));
      LiegroupElementRef rhs (space->elementRef (rhs_));
      vector_t error = output - space->neutral (); // log (rightHandSide)
      for (size_type k = 0; k < space->nv(); ++k) {
        if (comparison_[k] != Equality)
          assert (error[k] == 0);
      }
      rhs = space->neutral () + error; // exp (error)
    }

    vectorIn_t Implicit::rightHandSide () const
    {
      return rhs_;
    }

    vectorOut_t Implicit::rightHandSide ()
    {
      return rhs_;
    }

    void Implicit::rightHandSideFunction (const DifferentiableFunctionPtr_t& rhsF)
    {
      if (rhsF) {
        if (rhsF->inputSize() != 1 || rhsF->inputDerivativeSize() != 1) {
          HPP_THROW(std::invalid_argument, "Right hand side functions must take"
              "only one real value as input. Got " << rhsF->inputSize() <<
              " and " << rhsF->inputDerivativeSize());
        }
        if (*rhsF->outputSpace() != *function_->outputSpace()) {
          HPP_THROW(std::invalid_argument, "The right hand side function output"
              " space (" << *rhsF->outputSpace() << ") differs from the function"
              " output space (" << *function_->outputSpace() << ").");
        }
        // Check that the right hand side is non-constant on all axis.
        for (std::size_t i = 0; i < comparison_.size(); ++i) {
          if (comparison_[i] != constraints::Equality) {
            HPP_THROW(std::invalid_argument, "The comparison type of implicit "
                "constraint with right hand side function should be Equality "
                "on all dimensions.");
          }
        }
      }

      rhsFunction_ = rhsF;
    }

    vectorIn_t Implicit::rightHandSideAt (const value_type& s)
    {
      if (rhsFunction_) {
        vector_t S (1); S[0] = s;
        rhsFunction_->value (rhsFunction_->outputSpace()->elementRef(rhs_), S);
      }
      return rightHandSide();
    }

    size_type Implicit::parameterSize () const
    {
      return parameterSize_;
    }

    size_type Implicit::rightHandSideSize () const
    {
      return function_->outputSpace ()->nq ();
    }
    const ComparisonTypes_t& Implicit::comparisonType () const
    {
      return comparison_;
    }

    void Implicit::comparisonType (const ComparisonTypes_t& comp)
    {
      comparison_ = comp;
      parameterSize_ = computeParameterSize (comparison_);
    }

    Implicit::Implicit (const DifferentiableFunctionPtr_t& function,
                        ComparisonTypes_t comp) :
      comparison_ (comp), rhs_ (vector_t::Zero (function->outputSize ())),
      parameterSize_ (computeParameterSize (comparison_)),
      function_ (function)
    {
      if (comp.size () == 0) {
        // Argument was probably not provided, set to Equality
        comparison_ = ComparisonTypes_t (function->outputDerivativeSize (),
                                         Equality);
      }
      assert (function_->outputDerivativeSize () == (size_type)comparison_.size ());
    }

    Implicit::Implicit (const DifferentiableFunctionPtr_t& function,
                        ComparisonTypes_t comp, vectorIn_t rhs) :
      comparison_ (comp), rhs_ (rhs), function_ (function)
    {
      if (rhs_.size () == 0) {
        rhs_ = vector_t::Zero (function->outputSize ());
      }
      if (comp.size () == 0) {
        // Argument was probably not provided, set to Equality
        comparison_ = ComparisonTypes_t (function->outputDerivativeSize (),
                                         Equality);
      }
      parameterSize_ = computeParameterSize (comparison_);
      assert (function_->outputDerivativeSize () ==
              (size_type)comparison_.size ());
      assert (function->outputDerivativeSize () == (size_type) rhs_.size ());
    }

    Implicit::Implicit (const Implicit& other):
      comparison_ (other.comparison_), rhs_ (other.rhs_),
      parameterSize_ (other.parameterSize_), function_ (other.function_),
      rhsFunction_ (other.rhsFunction_)
    {
    }

    bool Implicit::isEqual (const Implicit& other, bool swapAndTest)
      const
    {
      if (comparison_ != other.comparison_) return false;
      if (rhs_.size() != other.rhs_.size()) return false;
      if (function_ != other.function_) return false;
      if (swapAndTest) return other.isEqual (*this, false);
      return true;
    }

    ImplicitPtr_t Implicit::create (
        const DifferentiableFunctionPtr_t& function)
    {
      ComparisonTypes_t comp (function->outputDerivativeSize(), constraints::EqualToZero);

      Implicit* ptr = new Implicit (function, comp);
      ImplicitPtr_t shPtr (ptr);
      ImplicitWkPtr_t wkPtr (shPtr);
      ptr->init (wkPtr);
      return shPtr;
    }

    ImplicitPtr_t Implicit::create (
        const DifferentiableFunctionPtr_t& function, ComparisonTypes_t comp)
    {
      Implicit* ptr = new Implicit (function, comp);
      ImplicitPtr_t shPtr (ptr);
      ImplicitWkPtr_t wkPtr (shPtr);
      ptr->init (wkPtr);
      return shPtr;
    }

    ImplicitPtr_t Implicit::create (
        const DifferentiableFunctionPtr_t& function,
        ComparisonTypes_t comp, vectorIn_t rhs)
    {
      Implicit* ptr = new Implicit (function, comp, rhs);
      ImplicitPtr_t shPtr (ptr);
      ImplicitWkPtr_t wkPtr (shPtr);
      ptr->init (wkPtr);
      return shPtr;
    }

    ImplicitPtr_t Implicit::createCopy (const ImplicitPtr_t& other)
    {
      Implicit* ptr = new Implicit (*other);
      ImplicitPtr_t shPtr (ptr);
      ImplicitWkPtr_t wkPtr (shPtr);
      ptr->init (wkPtr);
      return shPtr;
    }

    ImplicitPtr_t Implicit::copy () const
    {
      return createCopy (weak_.lock ());
    }

    void Implicit::rightHandSideFromConfig (ConfigurationIn_t config)
    {
      if (parameterSize () > 0) {
        LiegroupSpacePtr_t space (function_->outputSpace());
        LiegroupElement output (space);
        LiegroupElementRef rhs    (space->elementRef(rhs_));

        function_->value (output, config);
        // d.error is used here as an intermediate storage. The value
        // computed is not exactly the error
        vector_t error = output - rhs;
        for (size_type k = 0; k < space->nv(); ++k) {
          if (comparison_[k] != Equality)
            error[k] = 0;
        }
        rhs += error;
      }
    }
  } // namespace constraints
} // namespace hpp
