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

namespace hpp {
  namespace constraints {

    /// \addtogroup constraints
    /// \{

    /// Identity function
    /// \f$ q_{out} = q_{in} \f$
    ///
    /// \todo should we handle specifically this function is the solvers ?
    class HPP_CONSTRAINTS_DLLAPI Identity
      : public constraints::DifferentiableFunction
    {
      public:
        static IdentityPtr_t create (const LiegroupSpacePtr_t space, const std::string& name)
        {
          IdentityPtr_t ptr (new Identity (space, name));
          return ptr;
        }

        Identity (const LiegroupSpacePtr_t space, const std::string& name) :
          DifferentiableFunction (space->nq(), space->nv(), space, name) {}

      protected:
        void impl_compute (LiegroupElementRef y, vectorIn_t arg) const
        {
          y.vector() = arg;
        }

        void impl_jacobian (matrixOut_t J, vectorIn_t) const
        {
          J.setIdentity();
        }

        bool isEqual(const DifferentiableFunction& other) const {
          dynamic_cast<const Identity&>(other);
          if (!DifferentiableFunction::isEqual(other))
            return false;

          return true;
        }

      private:
        Identity() {}
        HPP_SERIALIZABLE();
    }; // class Identity

    /// Affine function
    /// \f$ f(q) = J * q + b \f$
    ///
    /// \todo should we handle specifically this function is the solvers ?
    class HPP_CONSTRAINTS_DLLAPI AffineFunction
      : public DifferentiableFunction
    {
    public:
      static AffineFunctionPtr_t create
        (const matrixIn_t& J, const std::string name = "LinearFunction")
      {
        return AffineFunctionPtr_t(new AffineFunction(J, name));
      }

      static AffineFunctionPtr_t create
        (const matrixIn_t& J, const vectorIn_t& b,
         const std::string name = "LinearFunction")
      {
        return AffineFunctionPtr_t(new AffineFunction(J, b, name));
      }
      protected:
        AffineFunction (const matrixIn_t& J,
            const std::string name = "LinearFunction")
          : DifferentiableFunction (J.cols(), J.cols(), LiegroupSpace::Rn
                                    (J.rows()), name),
          J_ (J), b_ (vector_t::Zero(J.rows()))
        {
          init();
        }

        AffineFunction (const matrixIn_t& J, const vectorIn_t& b,
            const std::string name = "LinearFunction")
          : DifferentiableFunction (J.cols(), J.cols(), LiegroupSpace::Rn
                                    (J.rows()), name),
          J_ (J), b_ (b)
        {
          init();
        }

        bool isEqual(const DifferentiableFunction& other) const {
          const AffineFunction& castother = dynamic_cast<const AffineFunction&>(other);
          if (!DifferentiableFunction::isEqual(other))
            return false;
          
          if (J_ != castother.J_)
            return false;
          if (b_ != castother.b_)
            return false;
          
          return true;
        }

      private:
        /// User implementation of function evaluation
        void impl_compute (LiegroupElementRef y, vectorIn_t x) const
        {
          y.vector ().noalias() = J_ * x + b_;
        }

        void impl_jacobian (matrixOut_t jacobian, vectorIn_t) const
        {
          jacobian = J_;
        }

        void init ()
        {
          assert(J_.rows() == b_.rows());
          activeParameters_ = (J_.array() != 0).colwise().any();
          activeDerivativeParameters_ = activeParameters_;
        }

        const matrix_t J_;
        const vector_t b_;

        AffineFunction() {}
        HPP_SERIALIZABLE();
    }; // class AffineFunction

    /// Constant function
    /// \f$ f(q) = C \f$
    ///
    /// \todo should we handle specifically this function is the solvers ?
    struct HPP_CONSTRAINTS_DLLAPI ConstantFunction
      : public DifferentiableFunction
    {
    public:
      static ConstantFunctionPtr_t create
        (const vector_t& constant, const size_type& sizeIn,
         const size_type& sizeInDer, const std::string name="ConstantFunction")
      {
        return ConstantFunctionPtr_t(new ConstantFunction(constant, sizeIn,
                                                          sizeInDer, name));
      }
      static ConstantFunctionPtr_t create
        (const LiegroupElement& element, const size_type& sizeIn,
         const size_type& sizeInDer, const std::string name="ConstantFunction")
      {
        return ConstantFunctionPtr_t(new ConstantFunction(element, sizeIn,
                                                          sizeInDer, name));
      }
    protected:
        ConstantFunction (const vector_t& constant,
                          const size_type& sizeIn,
                          const size_type& sizeInDer,
                          const std::string name = "ConstantFunction") :
          DifferentiableFunction (sizeIn, sizeInDer, LiegroupSpace::Rn
                                  (constant.rows()), name),
          c_ (constant, outputSpace())
        {}

        ConstantFunction (const LiegroupElement& element,
                          const size_type& sizeIn,
                          const size_type& sizeInDer,
                          const std::string name = "ConstantFunction") :
          DifferentiableFunction (sizeIn, sizeInDer, element.space(), name),
          c_ (element)
        {}

        /// User implementation of function evaluation
        void impl_compute (LiegroupElementRef r, vectorIn_t) const { r = c_; }

        void impl_jacobian (matrixOut_t J, vectorIn_t) const { J.setZero(); }

        bool isEqual(const DifferentiableFunction& other) const {
          const ConstantFunction& castother = dynamic_cast<const ConstantFunction&>(other);
          if (!DifferentiableFunction::isEqual(other))
            return false;
          
          if (c_.vector() == castother.c_.vector())
            return false;
          
          return true;
        }

        const LiegroupElement c_;

    private:
        ConstantFunction () {}
        HPP_SERIALIZABLE();
    }; // class ConstantFunction

    /// \}
  } // namespace constraints
} // namespace hpp


#endif // HPP_CONSTRAINTS_AFFINE_FUNCTION_HH
