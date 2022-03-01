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

#ifndef HPP_CONSTRAINTS_FUNCTION_OF_PARAMETER_SUBSET_HH
# define HPP_CONSTRAINTS_FUNCTION_OF_PARAMETER_SUBSET_HH

# include <hpp/constraints/differentiable-function.hh>

namespace hpp {
  namespace constraints {
    namespace function {
      /// Function depending on a subset of parameters
      ///
      /// This class implements a function over a configuration space,
      /// the output values of which only depend on a convex subset of
      /// parameters.
      ///
      /// \f{equation}
      /// f (q_1,\cdots,q_n) = g (q_{i},\cdots,q_{i+p-1})
      /// \f}
      /// where
      /// \li \f$n\f$ is the dimension of the input configuration space,
      /// \li \f$[i,i+p-1]\f$ is an interval included in \f$[1,n]\f$,
      /// \li \f$g\f$ is a differentiable mapping from \f$\mathbf{R}^p\f$
      ///     to the output space of \f$f\f$.
      class HPP_CONSTRAINTS_DLLAPI OfParameterSubset :
        public DifferentiableFunction
      {
      public:
        /// Create instance and return shared pointer
        /// \param g the mapping from the subset of parameters to the
        ///        output space,
        /// \param nArgs dimension \f$n\f$ of the input space representation,
        /// \param nDers dimension of the input tangent space,
        /// \param inArgs interval \f$[i,i+p-1]\f$ of configuration indices,
        /// \param inDers interval of velocity indices.
        static OfParameterSubsetPtr_t create
          (const DifferentiableFunctionPtr_t& g,
           const size_type& nArgs, const size_type& nDers,
           const segment_t& inArgs, const segment_t& inDers)
        {
          return OfParameterSubsetPtr_t
            (new OfParameterSubset (g, nArgs, nDers, inArgs, inDers));
        }

      protected:
        /// Constructor
        /// \param g the mapping from the subset of parameters to the
        ///        output space,
        /// \param nArgs dimension \f$n\f$ of the input space representation,
        /// \param nDers dimension of the input tangent space,
        /// \param inArgs interval \f$[i,i+p-1]\f$ of configuration indices,
        /// \param inDers interval of velocity indices.
        OfParameterSubset (const DifferentiableFunctionPtr_t& g,
                            const size_type& nArgs, const size_type& nDers,
                            const segment_t& inArgs, const segment_t& inDers);

        void impl_compute (LiegroupElementRef y, vectorIn_t arg) const
        {
          g_->value(y, arg.segment (sa_.first, sa_.second));
        }

        void impl_jacobian (matrixOut_t J, vectorIn_t arg) const
        {
          g_->jacobian(J.middleCols (sd_.first, sd_.second),
                           arg.segment (sa_.first, sa_.second));
        }

        bool isEqual(const DifferentiableFunction& other) const {
          const OfParameterSubset& castother = dynamic_cast<const OfParameterSubset&>(other);
          if (!DifferentiableFunction::isEqual(other))
            return false;
          
          if (g_ != castother.g_)
            return false;
          if (sa_ != castother.sa_)
            return false;
          if (sd_ != castother.sd_)
            return false;
          
          return true;
        }

        std::ostream& print (std::ostream& os) const;

        DifferentiableFunctionPtr_t g_;
        const segment_t sa_, sd_;
      }; // class OfParameterSubset
    } // namespace function
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_FUNCTION_OF_PARAMETER_SUBSET_HH
