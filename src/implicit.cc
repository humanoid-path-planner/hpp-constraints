// Copyright (c) 2015, LAAS-CNRS
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

#include "hpp/constraints/implicit.hh"

#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/weak_ptr.hpp>
#include <hpp/constraints/differentiable-function.hh>
#include <hpp/util/serialization.hh>
#include <pinocchio/serialization/eigen.hpp>

namespace hpp {
namespace constraints {
bool Implicit::operator==(const Implicit& other) const {
  bool res = isEqual(other, true);
  return res;
}
size_type computeParameterSize(const ComparisonTypes_t& comp) {
  size_type size = 0;
  for (std::size_t i = 0; i < comp.size(); ++i) {
    if (comp[i] == Equality) {
      ++size;
    }
  }
  return size;
}

void Implicit::rightHandSideFunction(const DifferentiableFunctionPtr_t& rhsF) {
  if (rhsF) {
    if (rhsF->inputSize() != 1 || rhsF->inputDerivativeSize() != 1) {
      HPP_THROW(std::invalid_argument,
                "Right hand side functions must take"
                "only one real value as input. Got "
                    << rhsF->inputSize() << " and "
                    << rhsF->inputDerivativeSize());
    }
    if (*rhsF->outputSpace() != *function_->outputSpace()) {
      HPP_THROW(std::invalid_argument,
                "The right hand side function output"
                " space ("
                    << *rhsF->outputSpace()
                    << ") differs from the function"
                       " output space ("
                    << *function_->outputSpace() << ").");
    }
    // Check that the right hand side is non-constant on all axis.
    for (std::size_t i = 0; i < comparison_.size(); ++i) {
      if (comparison_[i] != constraints::Equality) {
        HPP_THROW(std::invalid_argument,
                  "The comparison type of implicit "
                  "constraint with right hand side function should be Equality "
                  "on all dimensions.");
      }
    }
  }

  rhsFunction_ = rhsF;
}

vectorIn_t Implicit::rightHandSideAt(const value_type& s) {
  if (rhsFunction_) {
    vector_t S(1);
    S[0] = s;
    rhsFunction_->value(output_, S);
  }
  return output_.vector();
}

size_type Implicit::parameterSize() const { return parameterSize_; }

size_type Implicit::rightHandSideSize() const {
  return function_->outputSpace()->nq();
}
const ComparisonTypes_t& Implicit::comparisonType() const {
  return comparison_;
}

void Implicit::comparisonType(const ComparisonTypes_t& comp) {
  comparison_ = comp;
  parameterSize_ = computeParameterSize(comparison_);
  computeIndices();
}

void Implicit::setInactiveRowsToZero(vectorOut_t error) const {
  inactiveRows_.lview(error).setZero();
}

Implicit::Implicit(const DifferentiableFunctionPtr_t& function,
                   ComparisonTypes_t comp, std::vector<bool> mask)
    : comparison_(comp),
      rhs_(vector_t::Zero(function->outputSize())),
      parameterSize_(computeParameterSize(comparison_)),
      function_(function),
      mask_(mask),
      activeRows_(),
      inactiveRows_(),
      inequalityIndices_(),
      equalityIndices_(),
      output_(function->outputSpace()),
      logOutput_(function->outputSpace()->nv())

{
  // This constructor used to set comparison types to Equality if an
  // empty vector was given as input. Now you should provide the
  // comparison type at construction.
  assert(function_->outputDerivativeSize() == (size_type)comparison_.size());
  assert(function_->outputDerivativeSize() == (size_type)mask.size());
  computeActiveRows();
  computeIndices();
}

// Compute active rows
void Implicit::computeActiveRows() {
  segments_t inactiveRows;
  for (std::size_t i = 0; i < mask_.size(); ++i) {
    if (mask_[i])
      activeRows_.push_back(segment_t(i, 1));
    else
      inactiveRows.push_back(segment_t(i, 1));
  }
  Eigen::BlockIndex::shrink(activeRows_);
  Eigen::BlockIndex::shrink(inactiveRows);
  inactiveRows_ = Eigen::RowBlockIndices(inactiveRows);
}

void Implicit::computeIndices() {
  inequalityIndices_.clear();
  equalityIndices_.clearRows();
  for (std::size_t i = 0; i < comparison_.size(); ++i) {
    if ((comparison_[i] == Superior) || (comparison_[i] == Inferior)) {
      inequalityIndices_.push_back((size_type)i);
    } else if (comparison_[i] == Equality) {
      equalityIndices_.addRow((size_type)i, 1);
    }
  }
  equalityIndices_.updateRows<true, true, true>();
}

Implicit::Implicit(const Implicit& other)
    : comparison_(other.comparison_),
      rhs_(other.rhs_),
      parameterSize_(other.parameterSize_),
      function_(other.function_),
      rhsFunction_(other.rhsFunction_),
      mask_(other.mask_),
      activeRows_(other.activeRows_),
      inactiveRows_(other.inactiveRows_),
      inequalityIndices_(other.inequalityIndices_),
      equalityIndices_(other.equalityIndices_),
      output_(other.output_),
      logOutput_(other.logOutput_)

{}

bool Implicit::isEqual(const Implicit& other, bool swapAndTest) const {
  if (comparison_ != other.comparison_) return false;
  if (rhs_.size() != other.rhs_.size()) return false;
  if (function_ != other.function_ && *function_ != *other.function_)
    return false;
  if (swapAndTest) return other.isEqual(*this, false);
  return true;
}

ImplicitPtr_t Implicit::create(const DifferentiableFunctionPtr_t& function,
                               ComparisonTypes_t comp, std::vector<bool> mask) {
  if (mask.empty()) {
    mask = std::vector<bool>(function->outputSpace()->nv(), true);
  }
  Implicit* ptr = new Implicit(function, comp, mask);
  ImplicitPtr_t shPtr(ptr);
  ImplicitWkPtr_t wkPtr(shPtr);
  ptr->init(wkPtr);
  return shPtr;
}

ImplicitPtr_t Implicit::createCopy(const ImplicitPtr_t& other) {
  Implicit* ptr = new Implicit(*other);
  ImplicitPtr_t shPtr(ptr);
  ImplicitWkPtr_t wkPtr(shPtr);
  ptr->init(wkPtr);
  return shPtr;
}

ImplicitPtr_t Implicit::copy() const { return createCopy(weak_.lock()); }

void Implicit::rightHandSideFromConfig(ConfigurationIn_t config,
                                       LiegroupElementRef rhs) {
  assert(*rhs.space() == *function_->outputSpace());
  function_->value(output_, config);
  logOutput_.setZero();
  equalityIndices_.lview(logOutput_) = equalityIndices_.rview(log(output_));
  // rhs = exp(logOutput_)
  rhs = rhs.space()->exp(logOutput_);
}

bool Implicit::checkRightHandSide(LiegroupElementConstRef rhs) const {
  if (rhs.space() != function_->outputSpace()) return false;
  logOutput_ = log(rhs);
  for (std::size_t i = 0; i < comparison_.size(); ++i)
    if ((comparison_[i] != Equality) && (logOutput_[i] != 0)) return false;

  return true;
}

std::pair<JointConstPtr_t, JointConstPtr_t>
Implicit::doesConstrainRelPoseBetween(DeviceConstPtr_t robot) const {
  std::pair<JointConstPtr_t, JointConstPtr_t> joints =
      function_->dependsOnRelPoseBetween(robot);
  // check that the constraint fully constrains the relative pose
  if (function_->outputSpace()->nv() != 6 || !checkAllRowsActive()) {
    hppDout(
        info, "Constraint "
                  << function_->name()
                  << " does not fully constrain the relative pose of joints");
    return std::pair<JointConstPtr_t, JointConstPtr_t>(nullptr, nullptr);
  }
  return joints;
}

template <class Archive>
void Implicit::serialize(Archive& ar, const unsigned int version) {
  (void)version;
  ar& BOOST_SERIALIZATION_NVP(comparison_);
  ar& BOOST_SERIALIZATION_NVP(rhs_);
  if (!Archive::is_saving::value)
    parameterSize_ = computeParameterSize(comparison_);
  ar& BOOST_SERIALIZATION_NVP(function_);
  ar& BOOST_SERIALIZATION_NVP(rhsFunction_);
  ar& BOOST_SERIALIZATION_NVP(mask_);
  ar& BOOST_SERIALIZATION_NVP(weak_);
  if (!Archive::is_saving::value) {
    computeActiveRows();
    computeIndices();
  }
}

HPP_SERIALIZATION_IMPLEMENT(Implicit);
}  // namespace constraints
}  // namespace hpp
