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
#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include <hpp/model/center-of-mass-computation.hh>

namespace hpp {
  namespace constraints {
    namespace {
      static void convert (const constraints::vector3_t& v, model::vectorOut_t res)
      {
        res [0] = v[0]; res [1] = v[1]; res [2] = v[2];
      }

      static void convert (const constraints::matrix3_t& m, eigen::matrix3_t& res)
      {
        res (0,0) = m (0,0); res (0,1) = m (0,1); res (0,2) = m (0,2);
        res (1,0) = m (1,0); res (1,1) = m (1,1); res (1,2) = m (1,2);
        res (2,0) = m (2,0); res (2,1) = m (2,1); res (2,2) = m (2,2);
      }

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
        const JointPtr_t& jointRef, const vector3_t reference,
        std::vector <bool> mask)
    {
      CenterOfMassComputationPtr_t comc =
        CenterOfMassComputation::create (robot);
      comc->add (robot->rootJoint ());
      comc->computeMass ();
      return create (name, robot, comc,
                     jointL, jointR, jointRef, reference, mask);
    }

    ComBetweenFeetPtr_t ComBetweenFeet::create (
        const std::string& name, const DevicePtr_t& robot,
        const CenterOfMassComputationPtr_t& comc,
        const JointPtr_t& jointL, const JointPtr_t& jointR,
        const JointPtr_t& jointRef, const vector3_t reference,
        std::vector <bool> mask)
    {
      ComBetweenFeet* ptr = new ComBetweenFeet
        (name, robot, comc, jointL, jointR, jointRef, reference, mask);
      ComBetweenFeetPtr_t shPtr (ptr);
      return shPtr;
    }

    ComBetweenFeet::ComBetweenFeet (
        const std::string& name, const DevicePtr_t& robot,
        const CenterOfMassComputationPtr_t& comc,
        const JointPtr_t& jointL, const JointPtr_t& jointR,
        const JointPtr_t& jointRef, const vector3_t reference,
        std::vector <bool> mask) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
          size (mask), name),
      robot_ (robot), comc_ (comc), jointL_ (jointL), jointR_ (jointR),
      jointRef_ (jointRef), reference_ (reference), mask_ (mask),
      nominalCase_ (false), result_ (3), jacobian_ (3, robot->numberDof ())
    {
      cross_.setZero ();
      if (mask[0] && mask[1] && mask[2])
        nominalCase_ = true;
      jacobian_.setZero ();
    }

    void ComBetweenFeet::impl_compute (vectorOut_t result,
        ConfigurationIn_t argument)
      const throw ()
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      comc_->compute (Device::COM);
      const Transform3f& Ml = jointL_->currentTransformation ();
      const Transform3f& Mr = jointR_->currentTransformation ();
      const Transform3f& Mref = jointRef_->currentTransformation ();
      const vector3_t& x = comc_->com ();
      fcl::Matrix3f RT = Mref.getRotation (); RT.transpose ();

      const fcl::Vec3f& pl = Ml.getTranslation ();
      const fcl::Vec3f& pr = Mr.getTranslation ();
      const fcl::Vec3f& center = (pl + pr)*0.5;
      if (nominalCase_)
        convert (RT * (x - center) - reference_, result);
      else {
        convert (RT * (x - center) - reference_, result_);
        size_t index = 0;
        for (size_t i = 0; i < 3; ++i)
          if (mask_[i]) {
            result [index] = result_ [i];
            index++;
          }
      }
    }

    void ComBetweenFeet::impl_jacobian (matrixOut_t jacobian,
        ConfigurationIn_t arg) const throw ()
    {
      robot_->currentConfiguration (arg);
      robot_->computeForwardKinematics ();
      const ComJacobian_t& Jcom = comc_->jacobian ();
      const JointJacobian_t& Jl (jointL_->jacobian ());
      const JointJacobian_t& Jr (jointR_->jacobian ());
      const JointJacobian_t& Jref (jointRef_->jacobian ());
      const vector3_t& x = comc_->com ();

      const Transform3f& Ml = jointL_->currentTransformation ();
      const Transform3f& Mr = jointR_->currentTransformation ();
      const Transform3f& Mref = jointRef_->currentTransformation ();
      fcl::Matrix3f RT = Mref.getRotation (); RT.transpose ();
      eigen::matrix3_t eigenRT; convert (RT, eigenRT);
      const vector3_t& xr (Mr.getTranslation ());
      const vector3_t& xl (Ml.getTranslation ());
      const vector3_t c = (xl + xr) / 2;
      cross_ (0,1) = -x [2] + c [2]; cross_ (1,0) =  x [2] - c [2];
      cross_ (0,2) =  x [1] - c [1]; cross_ (2,0) = -x [1] + c [1];
      cross_ (1,2) = -x [0] + c [0]; cross_ (2,1) =  x [0] - c [0];

      if (nominalCase_) {
        jacobian.leftCols (Jref.cols ()) = eigenRT * (
            Jcom - 0.5 * ( Jl.topRows (3) + Jr.topRows (3) )
            - cross_ * Jref.bottomRows (3));
        jacobian.rightCols (jacobian.cols () - Jref.cols ()).setZero ();
      } else {
        jacobian_.leftCols (Jref.cols ()) = eigenRT * (
            Jcom - 0.5 * ( Jl.topRows (3) + Jr.topRows (3) )
            - cross_ * Jref.bottomRows (3));
        jacobian_.rightCols (jacobian.cols () - Jref.cols ()).setZero ();
        size_t index = 0;
        for (size_t i = 0; i < 3; ++i)
          if (mask_[i]) {
            jacobian.row (index) = jacobian_.row (i);
            index++;
          }
      }
    }
  } // namespace constraints
} // namespace hpp
