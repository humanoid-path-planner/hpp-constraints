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
#include <hpp/model/fcl-to-eigen.hh>
#include <hpp/model/joint.hh>
#include <hpp/constraints/position.hh>

namespace hpp {
  namespace constraints {
    static matrix3_t identity () { matrix3_t R; R.setIdentity (); return R;}
    matrix3_t Position::I3 = identity ();

    PositionPtr_t Position::create (const std::string& name,
                                    const DevicePtr_t& robot,
				    const JointPtr_t& joint,
				    const vector3_t& pointInLocalFrame,
				    const vector3_t& targetInGlobalFrame,
				    const matrix3_t& rotation,
				    std::vector <bool> mask)
    {
      Position* ptr = new Position
        (name, robot, joint, pointInLocalFrame, targetInGlobalFrame, rotation, mask);
      PositionPtr_t shPtr (ptr);
      return shPtr;
    }

    PositionPtr_t Position::create (const DevicePtr_t& robot,
				    const JointPtr_t& joint,
				    const vector3_t& pointInLocalFrame,
				    const vector3_t& targetInGlobalFrame,
				    const matrix3_t& rotation,
				    std::vector <bool> mask)
    {
      Position* ptr = new Position
        ("Position", robot, joint, pointInLocalFrame, targetInGlobalFrame, rotation, mask);
      PositionPtr_t shPtr (ptr);
      return shPtr;
    }

    Position::Position (const std::string& name, const DevicePtr_t& robot,
                        const JointPtr_t& joint, const vector3_t& pointInLocalFrame,
                        const vector3_t& targetInGlobalFrame, const matrix3_t& rotation,
                        std::vector <bool> mask) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
		mask [0] + mask [1] + mask [2],	name),
      robot_ (robot), joint_ (joint), pointInLocalFrame_ (pointInLocalFrame),
      targetInGlobalFrame_ (targetInGlobalFrame), SBT_ (), result_ (3)
    {
      if (rotation.isIdentity () && mask [0] && mask [1] && mask [2]) {
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

    void Position::impl_compute (vectorOut_t result,
				 ConfigurationIn_t argument)
      const throw ()
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      const Transform3f& M = joint_->currentTransformation ();
      p_ = M.transform (pointInLocalFrame_);
      model::toEigen (targetInGlobalFrame_ - p_, result_);
      if (!nominalCase_) result = SBT_*result_;
      else result = result_;
    }

    void Position::impl_jacobian (matrixOut_t jacobian,
				  ConfigurationIn_t arg) const throw ()
    {
      robot_->currentConfiguration (arg);
      robot_->computeForwardKinematics ();
      const Transform3f& M = joint_->currentTransformation ();
      p_ = M.getRotation () * pointInLocalFrame_;
      cross_ (0,1) = -p_[2]; cross_ (1,0) = p_[2];
      cross_ (0,2) = p_[1]; cross_ (2,0) = -p_[1];
      cross_ (1,2) = -p_[0]; cross_ (2,1) = p_[0];
      const JointJacobian_t& Jjoint (joint_->jacobian ());
      if (nominalCase_) {
	jacobian.leftCols (Jjoint.cols ()) =
	  cross_*Jjoint.bottomRows (3) - Jjoint.topRows (3);
	jacobian.rightCols (jacobian.cols () - Jjoint.cols ()).setZero ();
      }
      else {
	jacobian.leftCols (Jjoint.cols ()) =
	  SBT_*(cross_*Jjoint.bottomRows (3) - Jjoint.topRows (3));
	jacobian.rightCols (jacobian.cols () - Jjoint.cols ()).setZero ();
      }
    }
  } // namespace constraints
} // namespace hpp
