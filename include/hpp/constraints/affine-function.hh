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
        AffineFunction (const matrixIn_t& J,
            const std::string name = "LinearFunction")
          : DifferentiableFunction (J.cols(), J.cols(), LiegroupSpace::Rn
                                    (J.rows()), name),
          J_ (J), b_ (vector_t::Zero(J.rows())),
          aIdx_ (0, J_.cols()), JIdx_ (0, J_.cols()), C_JIdx_ (0,0),
          qshort_ (J.cols())
        {
          init();
        }

        AffineFunction (const matrixIn_t& J, const vectorIn_t& b,
            const std::string name = "LinearFunction")
          : DifferentiableFunction (J.cols(), J.cols(), LiegroupSpace::Rn
                                    (J.rows()), name),
          J_ (J), b_ (b),
          aIdx_ (0, J_.cols()), JIdx_ (0, J_.cols()), C_JIdx_ (0,0),
          qshort_ (J.cols())
        {
          init();
        }

        AffineFunction (const matrixIn_t& J, const vectorIn_t& b,
            Eigen::RowBlockIndices& argSelection,
            Eigen::ColBlockIndices& derSelection,
            Eigen::ColBlockIndices& CderSelection,
            const std::string name = "AffineFunction") :
          DifferentiableFunction (J.cols(), J.cols(), LiegroupSpace::Rn
                                  (J.rows()), name),
          J_ (J), b_ (b),
          aIdx_ (argSelection), JIdx_ (derSelection), C_JIdx_ (CderSelection),
          qshort_ (argSelection.nbIndices())
        {
          init();
        }

      private:
        /// User implementation of function evaluation
        void impl_compute (LiegroupElement& result, vectorIn_t argument) const
        {
          qshort_ = aIdx_.rview(argument);
          result.vector ().noalias () = J_ * qshort_;
          result.vector () += b_;
        }

        void impl_jacobian (matrixOut_t jacobian,
            vectorIn_t) const
        {
          JIdx_  .lview(jacobian) = J_;
          C_JIdx_.lview(jacobian).setZero();
        }

        void init ()
        {
          assert(J_.rows() == b_.rows());
          activeParameters_ = (J_.array() != 0).colwise().any();
          activeDerivativeParameters_ = activeParameters_;
        }

        const matrix_t J_;
        const vector_t b_;
        const Eigen::RowBlockIndices aIdx_;
        const Eigen::ColBlockIndices JIdx_, C_JIdx_;
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
                          const std::string name = "ConstantFunction") :
          DifferentiableFunction (sizeIn, sizeInDer, LiegroupSpace::Rn
                                  (constant.rows()), name),
          c_ (constant, LiegroupSpace::Rn (constant.rows()))
        {}

        /// User implementation of function evaluation
        void impl_compute (LiegroupElement& r, vectorIn_t) const { r = c_; }

        void impl_jacobian (matrixOut_t J, vectorIn_t) const { J.setZero(); }

        const LiegroupElement c_;
    }; // class ConstantFunction

    /// \}
  } // namespace constraints
} // namespace hpp


#endif // HPP_CONSTRAINTS_AFFINE_FUNCTION_HH
