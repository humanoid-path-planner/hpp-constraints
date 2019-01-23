// Copyright (c) 2017 - 2018, CNRS
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr), Florent Lamiraux
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

        std::ostream& print (std::ostream& os) const;

        DifferentiableFunctionPtr_t g_;
        const segment_t sa_, sd_;
      }; // class OfParameterSubset
    } // namespace function
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_FUNCTION_OF_PARAMETER_SUBSET_HH
