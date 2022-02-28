// Copyright (c) 2017 - 2018, CNRS
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr), Florent Lamiraux
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

#include <hpp/constraints/function/of-parameter-subset.hh>

#include <hpp/util/indent.hh>

namespace hpp {
  namespace constraints {
    namespace function {
      namespace {
        std::string toStr (const segment_t& s)
        {
          std::ostringstream os;
          os << "[ " << s.first << ", " << s.first + s.second << " ]";
          return os.str();
        }
      }

      OfParameterSubset::OfParameterSubset (const DifferentiableFunctionPtr_t& g,
                          const size_type& nArgs, const size_type& nDers,
                          const segment_t& inArgs, const segment_t& inDers) :
        DifferentiableFunction (nArgs, nDers, g->outputSpace(),
                                g->name() + " | " + toStr(inArgs)),
        g_ (g), sa_ (inArgs), sd_ (inDers)
      {
        activeParameters_.setConstant(false);
        activeParameters_.segment(sa_.first, sa_.second)
          = g_->activeParameters();

        activeDerivativeParameters_.setConstant(false);
        activeDerivativeParameters_.segment(sd_.first, sd_.second)
          = g_->activeDerivativeParameters();
      }

      std::ostream& OfParameterSubset::print (std::ostream& os) const
      {
        constraints::DifferentiableFunction::print(os);
        return os << incindent << iendl << *g_ << decindent;
      }
    } // namespace function
  } // namespace constraints
} // namespace hpp
