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
    void Implicit::rightHandSide (vectorIn_t rhs)
    {
      rhs_ = rhs;
    }

    vectorIn_t Implicit::rightHandSide () const
    {
      return rhs_;
    }

    vectorOut_t Implicit::nonConstRightHandSide ()
    {
      return rhs_;
    }

    size_type Implicit::rhsSize () const
    {
      return rhs_.size ();
    }

    const ComparisonTypes_t& Implicit::comparisonType () const
    {
      return comparison_;
    }

    bool Implicit::constantRightHandSide () const
    {
      for (std::size_t i = 0; i < comparison_.size(); ++i)
        if (comparison_[i] == constraints::Equality)
          return false;
      return true;
    }

    void Implicit::comparisonType (const ComparisonTypes_t& comp)
    {
      comparison_ = comp;
      if (constantRightHandSide())
        rhs_ = vector_t ();
      else
        rhs_ = vector_t::Zero (rhsRealSize_);
    }

    Implicit::Implicit (const DifferentiableFunctionPtr_t& function,
                        ComparisonTypes_t comp) :
      comparison_ (comp), rhs_ (vector_t::Zero (function->outputSize ())),
      rhsRealSize_ (rhs_.size()), function_ (function),
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
      // TODO currently, I think it is not supported to have an implicit constraint
      // with function_->outputSize() != function_->outputDerivativeSize()
      // because we apply the comparison type in the same manner to the value
      // coeffs and to the Jacobian rows.
      // EqualToZero should not cause any problem.
      // Equality is ill-defined coefficient wise (as it is implemented).
      // Superior and Inferior are (I think) ill-defined.
      // We should apply the comparison on `function_->value() - rhs`
      // and on the Jacobian rows, because they have the same dimension.
      assert (function_->outputSize () == function_->outputDerivativeSize ());
      assert (function_->outputSize () == comparison_.size ());
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
        rhs_.resize (function->outputDerivativeSize ());
      }
      if (constantRightHandSide ())
        rhs_ = vector_t ();
      /// TODO See TODO in the constructor above.
      assert (function_->outputSize () == function_->outputDerivativeSize ());
      assert (function_->outputSize () == comparison_.size ());
    }

    Implicit::Implicit (const Implicit& other):
      comparison_ (other.comparison_), rhs_ (other.rhs_),
      rhsRealSize_ (other.rhsRealSize_), function_ (other.function_),
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
      ComparisonTypes_t comp (function->outputSize(), constraints::EqualToZero);

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
      if (rhsSize () > 0) {
        assert (*(function_->outputSpace ()) ==
                *(LiegroupSpace::Rn (function_->outputSize ())));
        LiegroupElement value (function_->outputSpace ());
        function_->value (value, config);
        rightHandSide (value.vector ());
      }
    }
  } // namespace constraints
} // namespace hpp
