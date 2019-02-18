// Copyright (c) 2014, LAAS-CNRS
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

#ifndef HPP_CONSTRAINTS_STATIC_STABILITY_HH
# define HPP_CONSTRAINTS_STATIC_STABILITY_HH

# include <vector>

# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/config.hh>
# include <hpp/constraints/deprecated.hh>

# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/symbolic-calculus.hh>

namespace hpp {
  namespace constraints {

    /// \addtogroup constraints
    /// \{
    class HPP_CONSTRAINTS_DLLAPI StaticStability :
      public DifferentiableFunction {
      public:
        static const value_type G;
        static const Eigen::Matrix <value_type, 6, 1> Gravity;

        struct Contact_t {
          JointPtr_t joint1, joint2;
          vector3_t point1, point2;
          vector3_t normal1, normal2;
        };
        typedef std::vector <Contact_t> Contacts_t;

        /// Constructor
        /// \param robot the robot the constraints is applied to,
        /// \param com COM of the object in the joint frame.
        StaticStability (const std::string& name, const DevicePtr_t& robot,
            const Contacts_t& contacts,
            const CenterOfMassComputationPtr_t& com);

        static StaticStabilityPtr_t create (
            const std::string& name,
            const DevicePtr_t& robot,
            const Contacts_t& contacts,
            const CenterOfMassComputationPtr_t& com);

        static StaticStabilityPtr_t create (
            const DevicePtr_t& robot,
            const Contacts_t& contacts,
            const CenterOfMassComputationPtr_t& com);

        MatrixOfExpressions<>& phi () {
          return phi_;
        }

      private:
        void impl_compute (LiegroupElementRef result,
                           ConfigurationIn_t argument) const;

        void impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const;

        static void findBoundIndex (vectorIn_t u, vectorIn_t v, 
            value_type& lambdaMin, size_type* iMin,
            value_type& lambdaMax, size_type* iMax);

        /// Return false if uMinus.isZero(), i which case v also zero (not computed).
        bool computeUminusAndV (vectorIn_t u, vectorOut_t uMinus,
            vectorOut_t v) const;

        void computeVDot (const ConfigurationIn_t arg, vectorIn_t uMinus, vectorIn_t S,
            matrixIn_t uDot, matrixOut_t uMinusDot, matrixOut_t vDot) const;

        void computeLambdaDot (vectorIn_t u, vectorIn_t v, const std::size_t i0,
            matrixIn_t uDot, matrixIn_t vDot, vectorOut_t lambdaDot) const;

        DevicePtr_t robot_;
        Contacts_t contacts_;
        CenterOfMassComputationPtr_t com_;

        typedef MatrixOfExpressions<eigen::vector3_t, JacobianMatrix> MoE_t;

        mutable MoE_t phi_;
        mutable vector_t u_, uMinus_, v_;
        mutable matrix_t uDot_, uMinusDot_, vDot_;
        mutable vector_t lambdaDot_; 
    };
    /// \}
  } // namespace constraints
} // namespace hpp

#endif //  HPP_CONSTRAINTS_STATIC_STABILITY_HH
