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
# include <hpp/constraints/tools.hh>
# include <hpp/constraints/symbolic-calculus.hh>

namespace hpp {
  namespace constraints {

    /**
     *  Constraint on the relative position of the center of mass
     *
     *  The value of the function is defined as the position of the center
     *  of mass in the reference frame of a joint.
     *
     *  \f{eqnarray*}
     *  \mathbf{f}(\mathbf{q}) &=&
     *  \left(\begin{array}{c}
     *    ( x_{com} - x_{ref} ) \cdot u_z \\
     *    ( R^T (e \wedge u) ) \cdot u_z \\
     *    ( x_{com} - x_L ) \cdot (u)\\
     *    ( x_{com} - x_R ) \cdot (u)\\
     *  \end{array}\right)
     *  \f}
     *  where
     *  \li \f$\mathbf{x}_{com}\f$ is the position of the center of mass,
     *  \li \f$\mathbf{x_L}\f$ is the position of the left joint,
     *  \li \f$\mathbf{x_R}\f$ is the position of the right joint,
     *  \li \f$\mathbf{x}_{ref}\f$ is the desired position of the center of mass
     *      expressed in reference joint frame.
     *  \li \f$ u = x_R - x_L \f$
     *  \li \f$ e = x_{com} - (\frac{x_L + x_R}{2})\f$
    **/
    class HPP_CONSTRAINTS_DLLAPI ComBetweenFeet : public DifferentiableFunction
    {
      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        /// Return a shared pointer to a new instance
        static ComBetweenFeetPtr_t create (
            const std::string& name, const DevicePtr_t& robot,
            const JointPtr_t& jointLeft, const JointPtr_t& jointRight,
            const vector3_t   pointLeft, const vector3_t   pointRight,
            const JointPtr_t& jointReference, const vector3_t pointRef,
            std::vector <bool> mask = boost::assign::list_of (true)(true)(true)(true));

        /// Return a shared pointer to a new instance
        static ComBetweenFeetPtr_t create (
            const std::string& name, const DevicePtr_t& robot,
            const CenterOfMassComputationPtr_t& comc,
            const JointPtr_t& jointLeft, const JointPtr_t& jointRight,
            const vector3_t   pointLeft, const vector3_t   pointRight,
            const JointPtr_t& jointReference, const vector3_t pointRef,
            std::vector <bool> mask = boost::assign::list_of (true)(true)(true)(true));

        virtual ~ComBetweenFeet () throw () {}

        ComBetweenFeet (const std::string& name, const DevicePtr_t& robot,
            const CenterOfMassComputationPtr_t& comc,
            const JointPtr_t& jointLeft, const JointPtr_t& jointRight,
            const vector3_t   pointLeft, const vector3_t   pointRight,
            const JointPtr_t& jointReference, const vector3_t pointRef,
            std::vector <bool> mask);

      protected:
        /// Compute value of error
        ///
        /// \param argument configuration of the robot,
        /// \retval result error vector
        virtual void impl_compute (LiegroupElementRef result,
                                   ConfigurationIn_t argument) const throw ();

        virtual void impl_jacobian (matrixOut_t jacobian,
            ConfigurationIn_t arg) const throw ();
      private:
        DevicePtr_t robot_;
        mutable Traits<PointCom>::Ptr_t com_;
        Traits<PointInJoint>::Ptr_t left_, right_;
        eigen::vector3_t pointRef_;
        JointPtr_t jointRef_;
        typedef Difference < PointCom, PointInJoint > DiffPCPiJ;
        typedef Difference < PointInJoint, PointInJoint > DiffPiJPiJ;
        typedef Sum < PointInJoint, PointInJoint > SumPiJPiJ;
        typedef CrossProduct < Difference < PointCom,
                                            ScalarMultiply < SumPiJPiJ > >,
                               DiffPiJPiJ > ECrossU_t;
        mutable Traits<DiffPCPiJ>::Ptr_t xmxl_, xmxr_;
        mutable Traits<DiffPiJPiJ>::Ptr_t u_;
        mutable Traits<ECrossU_t>::Ptr_t ecrossu_;
        mutable Traits<RotationMultiply <ECrossU_t> >::Ptr_t expr_;
        mutable Traits<CalculusBaseAbstract<value_type, RowJacobianMatrix > >::Ptr_t xmxlDotu_, xmxrDotu_;
        std::vector <bool> mask_;
        mutable eigen::matrix3_t cross_;
    }; // class ComBetweenFeet
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_COM_BETWEEN_FEET_HH
