//
// Copyright (c) 2015 CNRS
// Authors: Joseph Mirabel
//
//

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

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
                              LiegroupSpace::Rn (size (mask)), name),
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

    void ComBetweenFeet::impl_compute (LiegroupElementRef result,
        ConfigurationIn_t argument)
      const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      size_t index = 0;
      if (mask_[0]) {
        com_->invalidate ();
        com_->computeValue (argument);
        result.vector () [index++] = (com_->value () - pointRef_)[2];
      }
      if (mask_[1]) {
        expr_->invalidate ();
        expr_->computeValue (argument);
        result.vector () [index++] = expr_->value ()[2];
      }
      if (mask_[2]) {
        xmxlDotu_->invalidate ();
        xmxlDotu_->computeValue (argument);
        result.vector () [index++] =   xmxlDotu_->value()[0];
      }
      if (mask_[3]) {
        xmxrDotu_->invalidate ();
        xmxrDotu_->computeValue (argument);
        result.vector () [index  ] =   xmxrDotu_->value()[0];
      }
    }

    void ComBetweenFeet::impl_jacobian (matrixOut_t jacobian,
        ConfigurationIn_t arg) const
    {
      robot_->currentConfiguration (arg);
      robot_->computeForwardKinematics ();
      size_t index = 0;
      if (mask_[0]) {
        com_->invalidate ();
        com_->computeJacobian (arg);
        jacobian.row (index++).leftCols (jointRef_->jacobian ().cols ())
          = com_->jacobian ().row (2);
      }
      if (mask_[1]) {
        expr_->invalidate ();
        expr_->computeJacobian (arg);
        jacobian.row (index++).leftCols (jointRef_->jacobian ().cols ())
          = expr_->jacobian ().row (2);
      }
      if (mask_[2]) {
        xmxlDotu_->invalidate ();
        xmxlDotu_->computeJacobian (arg);
        jacobian.row (index++).leftCols (jointRef_->jacobian ().cols ())
          = xmxlDotu_->jacobian ();
      }
      if (mask_[3]) {
        xmxrDotu_->invalidate ();
        xmxrDotu_->computeJacobian (arg);
        jacobian.row (index  ).leftCols (jointRef_->jacobian ().cols ())
          = xmxrDotu_->jacobian ();
      }
    }
  } // namespace _constraints
} // namespace hpp
