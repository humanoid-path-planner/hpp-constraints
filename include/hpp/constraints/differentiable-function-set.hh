//
// Copyright (c) 2016 CNRS
// Authors: Joseph Mirabel
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

#ifndef HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_SET_HH
#define HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_SET_HH

#include <hpp/constraints/differentiable-function.hh>
#include <hpp/constraints/fwd.hh>

namespace hpp {
namespace constraints {
/// \addtogroup constraints
/// \{

/// Set of differentiable functions
///
/// This class also handles selection of cols of the output matrix.
class HPP_CONSTRAINTS_DLLAPI DifferentiableFunctionSet
    : public DifferentiableFunction {
 public:
  typedef std::vector<DifferentiableFunctionPtr_t> Functions_t;

  /// Return a shared pointer to a new instance
  ///
  /// \param name the name of the constraints,
  static DifferentiableFunctionSetPtr_t create(const std::string& name) {
    return DifferentiableFunctionSetPtr_t(new DifferentiableFunctionSet(name));
  }

  virtual ~DifferentiableFunctionSet() {}

  /// \name Function stack management
  /// \{

  /// Get the stack of functions
  const Functions_t& functions() const { return functions_; }

  void add(const DifferentiableFunctionPtr_t& func) {
    if (functions_.empty()) {
      inputSize_ = func->inputSize();
      inputDerivativeSize_ = func->inputDerivativeSize();
      activeParameters_ = func->activeParameters();
      activeDerivativeParameters_ = func->activeDerivativeParameters();
    } else {
      assert(inputSize_ == func->inputSize());
      assert(inputDerivativeSize_ == func->inputDerivativeSize());

      activeParameters_ = activeParameters_ || func->activeParameters();
      activeDerivativeParameters_ =
          activeDerivativeParameters_ || func->activeDerivativeParameters();
    }
    functions_.push_back(func);
    result_.push_back(LiegroupElement(func->outputSpace()));
    *outputSpace_ *= func->outputSpace();
  }

  /// The output columns selection of other is not taken into account.
  void merge(const DifferentiableFunctionSetPtr_t& other) {
    const Functions_t& functions = other->functions();
    for (Functions_t::const_iterator _f = functions.begin();
         _f != functions.end(); ++_f)
      add(*_f);
  }

  /// \}

  std::ostream& print(std::ostream& os) const;

  /// Constructor
  ///
  /// \param name the name of the constraints,
  DifferentiableFunctionSet(const std::string& name)
      : DifferentiableFunction(0, 0, 0, name) {}

  DifferentiableFunctionSet() : DifferentiableFunction(0, 0, 0, "Stack") {}

 protected:
  void impl_compute(LiegroupElementRef result, ConfigurationIn_t arg) const {
    size_type row = 0;
    std::size_t i = 0;
    for (Functions_t::const_iterator _f = functions_.begin();
         _f != functions_.end(); ++_f) {
      const DifferentiableFunction& f = **_f;
      f.impl_compute(result_[i], arg);
      assert(hpp::pinocchio::checkNormalized(result_[i]));
      result.vector().segment(row, f.outputSize()) = result_[i].vector();
      row += f.outputSize();
      ++i;
    }
    assert(hpp::pinocchio::checkNormalized(result));
  }
  void impl_jacobian(matrixOut_t jacobian, ConfigurationIn_t arg) const {
    size_type row = 0;
    for (Functions_t::const_iterator _f = functions_.begin();
         _f != functions_.end(); ++_f) {
      const DifferentiableFunction& f = **_f;
      f.impl_jacobian(jacobian.middleRows(row, f.outputDerivativeSize()), arg);
      row += f.outputDerivativeSize();
    }
  }

  bool isEqual(const DifferentiableFunction& other) const {
    const DifferentiableFunctionSet& castother =
        dynamic_cast<const DifferentiableFunctionSet&>(other);
    if (!DifferentiableFunction::isEqual(other)) return false;

    if (functions_ != castother.functions_) return false;
    if (result_.size() != castother.result_.size()) return false;
    for (std::size_t i = 0; i < result_.size(); i++)
      if (result_[i] != castother.result_[i]) return false;

    return true;
  }

 private:
  Functions_t functions_;
  mutable std::vector<LiegroupElement> result_;
};  // class DifferentiableFunctionSet
/// \}
}  // namespace constraints
}  // namespace hpp

#endif  // HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_SET_HH
