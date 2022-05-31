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

#include <hpp/constraints/function/difference.hh>
#include <hpp/constraints/matrix-view.hh>

namespace hpp {
namespace constraints {
namespace function {
namespace {
std::string toStr(const segment_t& s) {
  std::ostringstream os;
  os << pretty_print(s);
  return os.str();
}
}  // namespace

Difference::Difference(const DifferentiableFunctionPtr_t& inner,
                       const size_type& nArgs, const size_type& nDers,
                       const segment_t& lInArgs, const segment_t& lInDers,
                       const segment_t& rInArgs, const segment_t& rInDers)
    : DifferentiableFunction(
          nArgs, nDers, LiegroupSpace::Rn(inner->outputSpace()->nv()),
          inner->name() + " | " + toStr(lInArgs) + " - " + toStr(rInArgs)),
      inner_(inner),
      lsa_(lInArgs),
      lsd_(lInDers),
      rsa_(rInArgs),
      rsd_(rInDers),
      l_(inner->outputSpace()),
      r_(inner->outputSpace()) {
  activeParameters_.setConstant(false);
  activeParameters_.segment(lsa_.first, lsa_.second) =
      inner_->activeParameters();
  activeParameters_.segment(rsa_.first, rsa_.second) =
      inner_->activeParameters();

  activeDerivativeParameters_.setConstant(false);
  activeDerivativeParameters_.segment(lsd_.first, lsd_.second) =
      inner_->activeDerivativeParameters();
  activeDerivativeParameters_.segment(rsd_.first, rsd_.second) =
      inner_->activeDerivativeParameters();
}

std::ostream& Difference::print(std::ostream& os) const {
  constraints::DifferentiableFunction::print(os);
  return os << incindent << iendl << *inner_ << decindent;
}
}  // namespace function
}  // namespace constraints
}  // namespace hpp
