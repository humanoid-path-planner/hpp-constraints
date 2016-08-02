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

#include <hpp/constraints/com-between-feet.hh>

#include <hpp/util/debug.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/center-of-mass-computation.hh>

namespace hpp {
  namespace constraints {
    namespace {
      static size_type size (std::vector<bool> mask)
      {
        size_type res = 0;
        for (std::vector<bool>::iterator it = mask.begin ();
            it != mask.end (); ++it) {
          if (*it) ++res;
        }
        return res;
      }
    } // namespace 

    ComBetweenFeetPtr_t ComBetweenFeet::create (
        const std::string& name, const DevicePtr_t& robot,
        const JointPtr_t& jointL, const JointPtr_t& jointR,
        const vector3_t   pointL, const vector3_t   pointR,
        const JointPtr_t& jointRef, const vector3_t pointRef,
        std::vector <bool> mask)
    {
      CenterOfMassComputationPtr_t comc =
        CenterOfMassComputation::create (robot);
      comc->add (robot->rootJoint ());
      return create (name, robot, comc, jointL, jointR,
          pointL, pointR, jointRef, pointRef, mask);
    }

    ComBetweenFeetPtr_t ComBetweenFeet::create (
        const std::string& name, const DevicePtr_t& robot,
        const CenterOfMassComputationPtr_t& comc,
        const JointPtr_t& jointL, const JointPtr_t& jointR,
        const vector3_t   pointL, const vector3_t   pointR,
        const JointPtr_t& jointRef, const vector3_t pointRef,
        std::vector <bool> mask)
    {
      ComBetweenFeet* ptr = new ComBetweenFeet
        (name, robot, comc, jointL, jointR, pointL, pointR, jointRef,
         pointRef, mask);
      ComBetweenFeetPtr_t shPtr (ptr);
      return shPtr;
    }

    ComBetweenFeet::ComBetweenFeet (
        const std::string& name, const DevicePtr_t& robot,
        const CenterOfMassComputationPtr_t& comc,
        const JointPtr_t& jointL, const JointPtr_t& jointR,
        const vector3_t   pointL, const vector3_t   pointR,
        const JointPtr_t& jointRef, const vector3_t pointRef,
        std::vector <bool> mask) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
          size (mask), name),
      robot_ (robot),
      com_ (PointCom::create (comc)),
      left_ (PointInJoint::create(jointL, pointL)),
      right_ (PointInJoint::create(jointR, pointR)),
      pointRef_ (),
      jointRef_ (jointRef),
      mask_ (mask)
    {
      cross_.setZero ();
      u_ = right_ - left_;
      xmxl_ = com_ - left_;
      xmxr_ = com_ - right_;
      ecrossu_ = (com_ - (0.5 * (left_ + right_)))^(u_);
      // expr_ = RotationMultiply <ECrossU_t> (jointRef_, ecrossu_, true);
      expr_ = JointTranspose (jointRef_) * ecrossu_;
      xmxlDotu_ = xmxl_ * u_;
      xmxrDotu_ = xmxr_ * u_;
      for (int i=0; i<3; i++) pointRef_[i] = pointRef[i];
    }

    void ComBetweenFeet::impl_compute (vectorOut_t result,
        ConfigurationIn_t argument)
      const throw ()
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      size_t index = 0;
      if (mask_[0]) {
        com_->invalidate ();
        com_->computeValue ();
        result[index++] = (com_->value () - pointRef_)[2];
      }
      if (mask_[1]) {
        expr_->invalidate ();
        expr_->computeValue ();
        result[index++] = expr_->value ()[2];
      }
      if (mask_[2]) {
        xmxlDotu_->invalidate ();
        xmxlDotu_->computeValue ();
        result[index++] =   xmxlDotu_->value();
      }
      if (mask_[3]) {
        xmxrDotu_->invalidate ();
        xmxrDotu_->computeValue ();
        result[index  ] =   xmxrDotu_->value();
      }
    }

    void ComBetweenFeet::impl_jacobian (matrixOut_t jacobian,
        ConfigurationIn_t arg) const throw ()
    {
      robot_->currentConfiguration (arg);
      robot_->computeForwardKinematics ();
      size_t index = 0;
      if (mask_[0]) {
        com_->invalidate ();
        com_->computeJacobian ();
        jacobian.row (index++).leftCols (jointRef_->jacobian ().cols ())
          = com_->jacobian ().row (2);
      }
      if (mask_[1]) {
        expr_->invalidate ();
        expr_->computeJacobian ();
        jacobian.row (index++).leftCols (jointRef_->jacobian ().cols ())
          = expr_->jacobian ().row (2);
      }
      if (mask_[2]) {
        xmxlDotu_->invalidate ();
        xmxlDotu_->computeJacobian ();
        jacobian.row (index++).leftCols (jointRef_->jacobian ().cols ())
          = xmxlDotu_->jacobian ();
      }
      if (mask_[3]) {
        xmxrDotu_->invalidate ();
        xmxrDotu_->computeJacobian ();
        jacobian.row (index  ).leftCols (jointRef_->jacobian ().cols ())
          = xmxrDotu_->jacobian ();
      }
    }
  } // namespace _constraints
} // namespace hpp
