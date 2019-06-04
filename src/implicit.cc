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
    size_type computeRightHandSideSize (const ComparisonTypes_t& comp)
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
    void Implicit::rightHandSide (vectorIn_t rhs)
    {
      rhs_ = rhs;
    }

    vectorIn_t Implicit::rightHandSide () const
    {
      return rhs_;
    }

    vectorOut_t Implicit::rightHandSide ()
    {
      return rhs_;
    }

    vectorOut_t Implicit::nonConstRightHandSide ()
    {
      return rhs_;
    }

    void Implicit::rightHandSideFunction (const DifferentiableFunctionPtr_t& rhsF)
    {
      if (rhsF) {
        assert (rhsF->inputSize() == 1);
        assert (rhsF->inputDerivativeSize() == 1);
        assert (*rhsF->outputSpace() == *function_->outputSpace());
        // Check that the right hand side is non-constant on all axis.
        for (std::size_t i = 0; i < comparison_.size(); ++i)
          assert (comparison_[i] == constraints::Equality);
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

    size_type Implicit::rhsSize () const
    {
      return rhs_.size ();
    }

    size_type Implicit::parameterSize () const
    {
      return rhsRealSize_;
    }

    size_type Implicit::rightHandSideSize () const
    {
      return function_->outputSpace ()->nq ();
    }
    const ComparisonTypes_t& Implicit::comparisonType () const
    {
      return comparison_;
    }

    bool Implicit::constantRightHandSide () const
    {
      for (std::size_t i = 0; i < comparison_.size(); ++i)
        if (comparison_[i] == constraints::Equality) {
          assert (parameterSize_ != 0);
          return false;
        }
      assert (parameterSize_ == 0);
      return true;
    }

    void Implicit::comparisonType (const ComparisonTypes_t& comp)
    {
      comparison_ = comp;
      if (constantRightHandSide())
        rhs_ = vector_t ();
      else
        rhs_ = vector_t::Zero (rhsRealSize_);
      rhsRealSize_ = computeRightHandSideSize (comparison_);
    }

    Implicit::Implicit (const DifferentiableFunctionPtr_t& function,
                        ComparisonTypes_t comp) :
      comparison_ (comp), rhs_ (vector_t::Zero (function->outputSize ())),
      rhsRealSize_ (computeRightHandSideSize (comparison_)),
      function_ (function),
      value_ (function->outputSize ()),
      jacobian_ (function->outputDerivativeSize (),
                 function->inputDerivativeSize ())
    {
      if (comp.size () == 0) {
        // Argument was probably not provided, set to Equality
        comparison_ = ComparisonTypes_t (function->outputDerivativeSize (),
                                         Equality);
      }
      if (constantRightHandSide ())
        rhs_ = vector_t ();
      assert (function_->outputDerivativeSize () == (size_type)comparison_.size ());
    }

    Implicit::Implicit (const DifferentiableFunctionPtr_t& function,
                        ComparisonTypes_t comp, vectorIn_t rhs) :
      comparison_ (comp), rhs_ (rhs), rhsRealSize_ (rhs.size()),
      function_ (function), value_ (function->outputSize ()),
      jacobian_ (matrix_t (function->outputSize (),
                           function->inputDerivativeSize ()))
    {
      if (comp.size () == 0) {
        // Argument was probably not provided, set to Equality
        comparison_ = ComparisonTypes_t (function->outputDerivativeSize (),
                                         Equality);
      }
      if (rhs.size () == 0) {
        rhs_.resize (function->outputSize ());
      }
      if (constantRightHandSide ())
        rhs_ = vector_t ();
      rhsRealSize_ = computeRightHandSideSize (comparison_);
      assert (function_->outputDerivativeSize () == (size_type)comparison_.size ());
    }

    Implicit::Implicit (const Implicit& other):
      comparison_ (other.comparison_), rhs_ (other.rhs_),
      rhsRealSize_ (other.rhsRealSize_), function_ (other.function_),
      rhsFunction_ (other.rhsFunction_),
      value_ (other.value_), jacobian_ (other.jacobian_)
    {
    }

    bool Implicit::isEqual (const Implicit& other, bool swapAndTest)
      const
    {
      if (comparison_ != other.comparison_) return false;
      if (rhs_.size() != other.rhs_.size() || rhs_ != other.rhs_) return false;
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

    // deprecated
    void Implicit::rightHandSideFromConfig (ConfigurationIn_t config)
    {
      if (parameterSize () > 0) {
        LiegroupElement value (function_->outputSpace ());
        function_->value (value, config);
        rightHandSide (value.vector ());
      }
    }
  } // namespace constraints
} // namespace hpp
