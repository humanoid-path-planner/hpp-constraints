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

#include <hpp/constraints/relative-com.hh>

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
            it != mask.end (); ++it)
          if (*it) ++res;
        return res;
      }
    } // namespace

    RelativeComPtr_t RelativeCom::create (const DevicePtr_t& robot,
					  const JointPtr_t& joint,
					  const vector3_t reference,
                                          std::vector <bool> mask)
    {
      CenterOfMassComputationPtr_t comc =
        CenterOfMassComputation::create (robot);
      comc->add (robot->rootJoint ());
      comc->computeMass ();
      return create (robot, comc, joint, reference, mask);
    }

    RelativeComPtr_t RelativeCom::create (
        const DevicePtr_t& robot,
        const CenterOfMassComputationPtr_t& comc,
        const JointPtr_t& joint, const vector3_t reference,
        std::vector <bool> mask)
    {
      RelativeCom* ptr = new RelativeCom (robot, comc, joint, reference, mask);
      RelativeComPtr_t shPtr (ptr);
      return shPtr;
    }

    RelativeCom::RelativeCom (const DevicePtr_t& robot,
        const CenterOfMassComputationPtr_t& comc,
        const JointPtr_t& joint, const vector3_t reference,
        std::vector <bool> mask) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
                               size (mask), "RelativeCom"),
      robot_ (robot), comc_ (comc), joint_ (joint), reference_ (reference), mask_ (mask),
      nominalCase_ (false), result_ (3), jacobian_ (3, robot->numberDof ())
    {
      cross_.setZero ();
      if (mask[0] && mask[1] && mask[2])
        nominalCase_ = true;
      jacobian_.setZero ();
    }

    void RelativeCom::impl_compute (vectorOut_t result,
				    ConfigurationIn_t argument)
      const throw ()
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      comc_->compute (Device::COM);
      const Transform3f& M = joint_->currentTransformation ();
      const vector3_t& x = comc_->com ();
      fcl::Matrix3f RT = M.getRotation (); RT.transpose ();
      const fcl::Vec3f& t = M.getTranslation ();

      if (nominalCase_)
        convert (RT * (x - t) - reference_, result);
      else {
        convert (RT * (x - t) - reference_, result_);
        size_t index = 0;
        for (size_t i = 0; i < 3; ++i)
          if (mask_[i]) {
            result [index] = result_ [i];
            index++;
          }
      }
    }

    void RelativeCom::impl_jacobian (matrixOut_t jacobian,
				     ConfigurationIn_t arg) const throw ()
    {
      robot_->currentConfiguration (arg);
      robot_->computeForwardKinematics ();
      comc_->compute (Device::ALL);
      const ComJacobian_t& Jcom = comc_->jacobian ();
      const JointJacobian_t& Jjoint (joint_->jacobian ());
      const Transform3f& M = joint_->currentTransformation ();
      fcl::Matrix3f RT (M.getRotation ()); RT.transpose ();
      const vector3_t& x = comc_->com ();
      const vector3_t& t (M.getTranslation ());
      cross_ (0,1) = -x [2] + t [2]; cross_ (1,0) = x [2] - t [2];
      cross_ (0,2) = x [1] - t [1]; cross_ (2,0) = -x [1] + t [1];
      cross_ (1,2) = -x [0] + t [0]; cross_ (2,1) = x [0] - t [0];
      eigen::matrix3_t eigenRT; convert (RT, eigenRT);
      if (nominalCase_) {
        jacobian.leftCols (Jjoint.cols ()) =
          eigenRT * (Jcom + cross_ * Jjoint.bottomRows (3) - Jjoint.topRows (3));
        jacobian.rightCols (jacobian.cols () - Jjoint.cols ()).setZero ();
      } else {
        jacobian_.leftCols (Jjoint.cols ()) =
          eigenRT * (Jcom + cross_ * Jjoint.bottomRows (3) - Jjoint.topRows (3));
        jacobian_.rightCols (jacobian.cols () - Jjoint.cols ()).setZero ();
        size_t index = 0;
        for (size_t i = 0; i < 3; ++i)
          if (mask_[i]) {
            jacobian.row (index) = jacobian_.row (i);
            index++;
          }
      }
      hppDout (info, "Jcom = " << std::endl << Jcom);
      hppDout (info, "Jw = " << std::endl << Jjoint.bottomRows (3));
      hppDout (info, "Jv = " << std::endl << Jjoint.topRows (3));
    }

  } // namespace constraints
} // namespace hpp
