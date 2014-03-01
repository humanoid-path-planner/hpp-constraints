//
// Copyright (c) 2014 CNRS
// Authors: Florent Lamiraux
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

#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include <hpp/constraints/relative-position.hh>

namespace hpp {
  namespace constraints {
    extern void applyRigidBodyMotion (const matrix4d& M, const vector3d& v,
				      vector3d& p);
    extern void applyRotation (const matrix4d& M, const vector3d& v,
			       vector3d& p);

    static vector3d zero3d (0, 0, 0);
    RelativePositionPtr_t RelativePosition::create
    (const DevicePtr_t& robot, const JointPtr_t&  joint1,
     const JointPtr_t& joint2, const vector3d& point1,
     const vector3d& point2)
    {
      RelativePosition* ptr = new RelativePosition
	(robot, joint1, joint2, point1, point2);
      RelativePositionShPtr shPtr (ptr);
      return shPtr;
    }

    static void cross (const vector3d& v, matrix3d& m)
    {
      m (0,1) = -v [2]; m (1,0) = v [2];
      m (0,2) = v [1]; m (2,0) = -v [1];
      m (1,2) = -v [0]; m (2,1) = v [0];
    }

    RelativePosition::RelativePosition (const DevicePtr_t& robot,
					const JointPtr_t& joint1,
					const JointPtr_t& joint2,
					const vector3d& point1,
					const vector3d& point2) :
      parent_t (robot->numberDof (), 3, "RelativePosition"),
      robot_ (robot), joint1_ (joint1), joint2_ (joint2),
      point1_ (point1), point2_ (point2), Jjoint1_ (6, robot->numberDof ()),
      Jjoint2_ (6, robot->numberDof ())
    {
      cross1_.setZero ();
      cross2_.setZero ();
    }

     void RelativePosition::impl_compute
     (result_t& result, const argument_t& argument) const throw ()
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      const matrix4d& M1 = joint1_->currentTransformation ();
      const matrix4d& M2 = joint2_->currentTransformation ();
      applyRigidBodyMotion (M1, point1_, global1_);
      applyRigidBodyMotion (M2, point2_, global2_);
      result = global1_ - global2_;
    }

  void RelativePosition::impl_jacobian (jacobian_t &jacobian,
					const argument_t &arg) const throw ()
    {
      robot_->currentConfiguration (arg);
      robot_->computeForwardKinematics ();
      const matrix4d& M1 = joint1_->currentTransformation ();
      const matrix4d& M2 = joint2_->currentTransformation ();
      robot_->getJacobian (*(robot_->rootJoint ()), *(joint1_->dynamic ()),
			   zero3d, Jjoint1_);
      robot_->getJacobian (*(robot_->rootJoint ()), *(joint2_->dynamic ()),
			   zero3d, Jjoint2_);
      applyRotation (M1, point1_, R1x1_);
      // -[R1 (q) x1]x
      cross (-R1x1_, cross1_);
      applyRotation (M2, point2_, R2x2_);
      // [R2 (q) x2]x
      cross (R2x2_, cross2_);
      jacobian = cross1_ * Jjoint1_.bottomRows (3)
	+ Jjoint1_.topRows (3) + cross2_ * Jjoint2_.bottomRows (3)
	- Jjoint2_.topRows (3);
    }

    void RelativePosition::impl_gradient (gradient_t &gradient,
					  const argument_t &argument,
					  size_type functionId) const throw ()
    {
      matrix_t jacobian (outputSize (), inputSize ());
      impl_jacobian (jacobian, argument);
      gradient = jacobian.row (functionId);
    }

  } // namespace constraints
} // namespace hpp
