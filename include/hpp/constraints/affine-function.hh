// Copyright (c) 2017, Joseph Mirabel
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

#ifndef HPP_CONSTRAINTS_AFFINE_FUNCTION_HH
# define HPP_CONSTRAINTS_AFFINE_FUNCTION_HH

# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/config.hh>

# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/matrix-view.hh>

namespace hpp {
  namespace constraints {
    /// \addtogroup constraints
    /// \{

    /// Affine function
    /// \f$ f(q) = J * q + b \f$
    class HPP_CONSTRAINTS_DLLAPI AffineFunction
      : public DifferentiableFunction
    {
      public:
        AffineFunction (const matrix_t& J, const vector_t& b,
            Eigen::RowBlockIndexes& argSelection,
            Eigen::ColBlockIndexes& derSelection,
            Eigen::ColBlockIndexes& CderSelection,
            const std::string name = "AffineFunction")
          : DifferentiableFunction (J.cols(), J.cols(), J.rows(), J.rows(), name),
          J_ (J), b_ (b),
          aIdx_ (argSelection), JIdx_ (derSelection), C_JIdx_ (CderSelection),
          qshort_ (argSelection.nbIndexes())
          {}

      private:
        /// User implementation of function evaluation
        void impl_compute (vectorOut_t result,
            vectorIn_t argument) const
        {
          aIdx_.view(argument).writeTo (qshort_);
          result.noalias() = J_ * aIdx_;
          result += b_;
        }

        void impl_jacobian (matrixOut_t jacobian,
            vectorIn_t) const
        {
          JIdx_  .view(jacobian) = J_;
          C_JIdx_.view(jacobian).setZero();
        }

        const matrix_t J_;
        const vector_t b_;
        const Eigen::RowBlockIndexes aIdx_;
        const Eigen::ColBlockIndexes JIdx_, C_JIdx_;
        mutable vector_t qshort_;
    }; // class AffineFunction

    /// Constant function
    /// \f$ f(q) = C \f$
    struct HPP_CONSTRAINTS_DLLAPI ConstantFunction
      : public DifferentiableFunction
    {
        ConstantFunction (const vector_t& constant,
                          const size_type& sizeIn, 
                          const size_type& sizeInDer,
                          const std::string name = "ConstantFunction")
          : DifferentiableFunction (sizeIn, sizeInDer, constant.rows(), constant.rows(), name),
          c_ (constant)
        {}

        /// User implementation of function evaluation
        void impl_compute (vectorOut_t r, vectorIn_t) const { r = c_; }

        void impl_jacobian (matrixOut_t J, vectorIn_t) const { J.setZero(); }

        const vector_t c_;
    }; // class ConstantFunction

    /// \}
  } // namespace constraints
} // namespace hpp


#endif // HPP_CONSTRAINTS_AFFINE_FUNCTION_HH
