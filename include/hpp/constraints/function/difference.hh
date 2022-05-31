// Copyright (c) 2018, Joseph Mirabel
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

#ifndef HPP_CONSTRAINTS_FUNCTION_DIFFERENCE_HH
#define HPP_CONSTRAINTS_FUNCTION_DIFFERENCE_HH

#include <hpp/constraints/config.hh>
#include <hpp/constraints/differentiable-function.hh>
#include <hpp/constraints/fwd.hh>

namespace hpp {
namespace constraints {
namespace function {
class Difference;
typedef shared_ptr<Difference> DifferencePtr_t;

/// Compute the difference between the value of the function in two points.
/// i.e.: \f$ f (q_0, ... , q_n) = f_{inner} (q_{left}) - f_{inner} (q_{right})
/// \f$
class HPP_CONSTRAINTS_LOCAL Difference
    : public constraints::DifferentiableFunction {
 public:
  typedef shared_ptr<Difference> Ptr_t;
  Difference(const DifferentiableFunctionPtr_t& inner, const size_type& nArgs,
             const size_type& nDers, const segment_t& lInArgs,
             const segment_t& lInDers, const segment_t& rInArgs,
             const segment_t& rInDers);

 protected:
  void impl_compute(LiegroupElementRef y, vectorIn_t arg) const {
    inner_->value(l_, arg.segment(lsa_.first, lsa_.second));
    inner_->value(r_, arg.segment(rsa_.first, rsa_.second));
    y.vector() = l_ - r_;
  }

  void impl_jacobian(matrixOut_t J, vectorIn_t arg) const {
    inner_->jacobian(J.middleCols(lsd_.first, lsd_.second),
                     arg.segment(lsa_.first, lsa_.second));
    inner_->jacobian(J.middleCols(rsd_.first, rsd_.second),
                     arg.segment(rsa_.first, rsa_.second));
    J.middleCols(rsd_.first, rsd_.second) *= -1;
  }

  bool isEqual(const DifferentiableFunction& other) const {
    const Difference& castother = dynamic_cast<const Difference&>(other);
    if (!DifferentiableFunction::isEqual(other)) return false;

    if (inner_ != castother.inner_) return false;
    if (lsa_ != castother.lsa_) return false;
    if (lsd_ != castother.lsd_) return false;
    if (rsa_ != castother.rsa_) return false;
    if (rsd_ != castother.rsd_) return false;

    return true;
  }

  std::ostream& print(std::ostream& os) const;

  DifferentiableFunctionPtr_t inner_;
  const segment_t lsa_, lsd_;
  const segment_t rsa_, rsd_;

  mutable LiegroupElement l_, r_;
};  // class Difference
}  // namespace function
}  // namespace constraints
}  // namespace hpp

#endif  // HPP_CONSTRAINTS_FUNCTION_DIFFERENCE_HH
