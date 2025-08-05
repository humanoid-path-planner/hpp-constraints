//
// Copyright (c) 2025 CNRS
// Authors: Florent Lamiraux
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

#ifndef HPP_CONSTRAINTS__EXTRACTED_FUNCTION_HH
#define HPP_CONSTRAINTS__EXTRACTED_FUNCTION_HH

#include <hpp/constraints/differentiable-function.hh>

namespace hpp {
namespace constraints {

class ExtractedFunction : public DifferentiableFunction {
public:
  typedef shared_ptr<ExtractedFunction> Ptr_t;
  typedef weak_ptr<ExtractedFunction> WkPtr_t;
  static Ptr_t create(const DifferentiableFunctionPtr_t& original, interval_t paramRange) {
    ExtractedFunction* ptr(new ExtractedFunction(original, paramRange));
    Ptr_t shPtr(ptr);
    ptr->weak_ = shPtr;
    return shPtr;
  }
  virtual ~ExtractedFunction() {}
  virtual void impl_compute(LiegroupElementRef result, vectorIn_t argument) const {
    assert(argument.rows() == 1);
    assert(argument.cols() == 1);
    Eigen::Matrix<value_type, 1, 1> param;
    if (reversed_) {
      param[0] = t12_ - argument[0];
      return original_->value(result, param);
    } else {
      param[0] = t12_ + argument[0];
      return original_->value(result, param);
    }
  }

  virtual void impl_jacobian(matrixOut_t jacobian, vectorIn_t arg) const {
    assert(arg.rows() == 1);
    assert(arg.cols() == 1);
    Eigen::Matrix<value_type, 1, 1> param;
    if (reversed_) {
      param[0] = t12_ - arg[0];
      original_->jacobian(jacobian, param);
      jacobian *= -1;
    } else {
      param[0] = arg[0] - t12_;
      return original_->jacobian(jacobian, param);
    }
  }

private:
  ExtractedFunction(const DifferentiableFunctionPtr_t& original, interval_t paramRange) :
    DifferentiableFunction(original->inputSize(), original->inputDerivativeSize(),
			   original->outputSpace(),
			   std::string("extracted from ") + original->name()),
    reversed_(paramRange.first > paramRange.second),
    t12_(reversed_ ? paramRange.second : paramRange.first), original_(original) {
    if (inputSize_ != 1) {
      std::ostringstream os;
      os << "hpp::constraints::DifferentiableFunction::extract: input size (=" << inputSize_
	 << ") should be equal to 1.";
      throw std::logic_error(os.str().c_str());
    }
    if (inputDerivativeSize_ != 1) {
      std::ostringstream os;
      os << "hpp::constraints::DifferentiableFunction::extract: input derivative size (="
	 << inputDerivativeSize_ << ") should be equal to 1.";
      throw std::logic_error(os.str().c_str());
    }
  }
  bool reversed_;
  value_type t12_;
  DifferentiableFunctionPtr_t original_;
  WkPtr_t weak_;
}; // class ExtractedFunction

} // namespace constraints
} // namespace hpp
#endif
