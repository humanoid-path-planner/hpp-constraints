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

#include <vector>
#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include <hpp/constraints/position.hh>

namespace hpp {
  namespace constraints {

    static vector3d zero3d (0, 0, 0);
    void applyRigidBodyMotion (const matrix4d& M, const vector3d& v,
				      vector3d& p)
    {
      p =  M.topLeftCorner <3,3> () * v + M.topRightCorner <3,1> ();
    }

    void applyRotation (const matrix4d& M, const vector3d& v,
			vector3d& p)
    {
      p =  M.topLeftCorner <3,3> () * v;
    }

    PositionPtr_t Position::create (const DevicePtr_t& robot,
				    const JointPtr_t& joint,
				    const vector3d& pointInLocalFrame,
				    const vector3d& targetInGlobalFrame,
				    const matrix3d& rotation,
				    std::vector <bool> mask)
    {
      Position* ptr = new Position (robot, joint, pointInLocalFrame,
				    targetInGlobalFrame, rotation, mask);
      PositionPtr_t shPtr (ptr);
      return shPtr;
    }

    Position::Position (const DevicePtr_t& robot, const JointPtr_t& joint,
			const vector3d& pointInLocalFrame,
			const vector3d& targetInGlobalFrame,
			const matrix3d& rotation,
			std::vector <bool> mask) :
      parent_t (robot->numberDof (), mask [0] + mask [1] + mask [2],
		"Position"),
      robot_ (robot), joint_ (joint), pointInLocalFrame_ (pointInLocalFrame),
      targetInGlobalFrame_ (targetInGlobalFrame), SBT_ (),
      Jjoint_ (6, robot->numberDof ())
    {
      if (rotation == Eigen::Matrix<double, 3, 3>::Identity() &&
	  mask [0] && mask [1] && mask [2]) {
	nominalCase_ = true;
      } else {
	nominalCase_ = false;
	size_type nbRows = mask [0] + mask [1] + mask [2];
	matrix_t selection (nbRows, 3); selection.setZero ();
	size_type row = 0;
	for (std::size_t i=0; i<3; ++i) {
	  if (mask [i]) {
	    selection (row, i) = 1;
	    ++row;
	  }
	}
	SBT_ = selection * rotation;
      }
      cross_.setZero ();
    }

    void Position::impl_compute (result_t& result,
				 const argument_t& argument) const throw ()
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      const matrix4d& M = joint_->currentTransformation ();
      applyRigidBodyMotion (M, pointInLocalFrame_, p_);
      result = targetInGlobalFrame_ - p_;
      if (!nominalCase_) result = SBT_*result;
    }

    void Position::impl_jacobian (jacobian_t &jacobian,
				  const argument_t &arg) const throw ()
    {
      robot_->currentConfiguration (arg);
      robot_->computeForwardKinematics ();
      const matrix4d& M = joint_->currentTransformation ();
      applyRotation (M, pointInLocalFrame_, p_);
      cross_ (0,1) = -p_[2]; cross_ (1,0) = p_[2];
      cross_ (0,2) = p_[1]; cross_ (0,2) = -p_[1];
      cross_ (1,2) = -p_[0]; cross_ (2,1) = p_[0];
      robot_->getJacobian (*(robot_->rootJoint ()), *(joint_->dynamic ()),
			   zero3d, Jjoint_);
      if (nominalCase_)
	jacobian = cross_*Jjoint_.bottomRows (3) - Jjoint_.topRows (3);
      else
	jacobian = SBT_*(cross_*Jjoint_.bottomRows (3) - Jjoint_.topRows (3));

    }
    void Position::impl_gradient (gradient_t &gradient,
				  const argument_t &argument,
				  size_type functionId) const throw ()
    {
      matrix_t jacobian (outputSize (), inputSize ());
      impl_jacobian (jacobian, argument);
      gradient = jacobian.row (functionId);
    }

  } // namespace constraints
} // namespace hpp
