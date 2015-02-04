//
// Copyright (c) 2015 CNRS
// Authors: Joseph Mirabel
//
//
// This file is part of hpp-constraints.
// hpp-constraints is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-constraints is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-constraints. If not, see
// <http://www.gnu.org/licenses/>.

#ifndef HPP_CONSTRAINTS_COM_BETWEEN_FEET_HH
# define HPP_CONSTRAINTS_COM_BETWEEN_FEET_HH

# include <boost/assign/list_of.hpp>
# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/config.hh>
# include <hpp/constraints/fwd.hh>

namespace hpp {
  namespace constraints {
    /// Constraint on the relative position of the center of mass
    ///
    /// The value of the function is defined as the position of the center
    /// of mass in the reference frame of a joint.
    /**
     *  \f{eqnarray*}
     *  \mathbf{f}(\mathbf{q}) &=&
     *  R^T \left(\mathbf{x} - \frac{x_L + x_R}{2}\right) - \mathbf{x}^{*}\\
     *  \mathbf{\dot{f}} &=& R^T \left(
     *  J_{com} - \frac{J_L^{\mathbf{v}} + J_R^{\mathbf{v}}}{2}
     *  + [\mathbf{x}-\frac{x_L + x_R}{2}]_{\times}J_{ref}^{\omega}
     *  \right) \mathbf{\dot{q}}
     *  \f}
    **/
    /// where
    /// \li \f$
    ///     \left(\begin{array}{cc} R & \mathbf{t} \\ 0 & 1\end{array}\right)
    ///     \f$
    /// is the position of the reference joint,
    /// \li \f$\mathbf{x}\f$ is the position of the center of mass,
    /// \li \f$\mathbf{x_L}\f$ is the position of the left joint,
    /// \li \f$\mathbf{x_R}\f$ is the position of the right joint,
    /// \li \f$\mathbf{x}^{*}\f$ is the desired position of the center of mass
    ///     expressed in reference joint frame.
    class HPP_CONSTRAINTS_DLLAPI ComBetweenFeet : public DifferentiableFunction
    {
      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        /// Return a shared pointer to a new instance
        static ComBetweenFeetPtr_t create (
            const std::string& name, const DevicePtr_t& robot,
            const JointPtr_t& jointLeft, const JointPtr_t& jointRight,
            const JointPtr_t& jointReference, const vector3_t reference,
            std::vector <bool> mask = boost::assign::list_of (true)(true)(true));

        /// Return a shared pointer to a new instance
        static ComBetweenFeetPtr_t create (
            const std::string& name, const DevicePtr_t& robot,
            const CenterOfMassComputationPtr_t& comc,
            const JointPtr_t& jointLeft, const JointPtr_t& jointRight,
            const JointPtr_t& jointReference, const vector3_t reference,
            std::vector <bool> mask = boost::assign::list_of (true)(true)(true));

        virtual ~ComBetweenFeet () throw () {}

        ComBetweenFeet (const std::string& name, const DevicePtr_t& robot,
            const CenterOfMassComputationPtr_t& comc,
            const JointPtr_t& jointLeft, const JointPtr_t& jointRight,
            const JointPtr_t& jointReference, const vector3_t reference,
            std::vector <bool> mask);

      protected:
        /// Compute value of error
        ///
        /// \param argument configuration of the robot,
        /// \retval result error vector
        virtual void impl_compute (vectorOut_t result,
            ConfigurationIn_t argument) const throw ();

        virtual void impl_jacobian (matrixOut_t jacobian,
            ConfigurationIn_t arg) const throw ();
      private:
        DevicePtr_t robot_;
        CenterOfMassComputationPtr_t comc_;
        JointPtr_t jointL_, jointR_, jointRef_;
        vector3_t reference_;
        std::vector <bool> mask_;
        bool nominalCase_;
        mutable eigen::matrix3_t cross_;
        mutable vector_t result_;
        mutable ComJacobian_t jacobian_;
    }; // class ComBetweenFeet
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_COM_BETWEEN_FEET_HH
