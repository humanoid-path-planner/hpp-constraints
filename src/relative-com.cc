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
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/center-of-mass-computation.hh>
#include <hpp/pinocchio/liegroup-element.hh>

#include <hpp/constraints/macros.hh>

namespace hpp {
  namespace constraints {

    namespace {
      static size_type size (std::vector<bool> mask)
      {
        size_type res = 0;
        for (std::vector<bool>::iterator it = mask.begin ();
            it != mask.end (); ++it)
          if (*it) ++res;
        return res;
      }
    } // namespace

    RelativeComPtr_t RelativeCom::create (const std::string& name,
                                          const DevicePtr_t& robot,
					  const JointPtr_t& joint,
					  const vector3_t reference,
                                          std::vector <bool> mask)
    {
      CenterOfMassComputationPtr_t comc =
        CenterOfMassComputation::create (robot);
      comc->add (robot->rootJoint ());
      return create (name, robot, comc, joint, reference, mask);
    }

    RelativeComPtr_t RelativeCom::create (
        const DevicePtr_t& robot,
        const CenterOfMassComputationPtr_t& comc,
        const JointPtr_t& joint, const vector3_t reference,
        std::vector <bool> mask)
    {
      return create ("RelativeCom", robot, comc, joint, reference, mask);
    }

    RelativeComPtr_t RelativeCom::create (
        const std::string& name,
        const DevicePtr_t& robot,
        const CenterOfMassComputationPtr_t& comc,
        const JointPtr_t& joint, const vector3_t reference,
        std::vector <bool> mask)
    {
      RelativeCom* ptr = new RelativeCom (robot, comc, joint, reference, mask, name);
      RelativeComPtr_t shPtr (ptr);
      return shPtr;
    }

    RelativeCom::RelativeCom (const DevicePtr_t& robot,
        const CenterOfMassComputationPtr_t& comc,
        const JointPtr_t& joint, const vector3_t reference,
        std::vector <bool> mask,
        const std::string& name) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
                              LiegroupSpace::Rn (size (mask)), name),
      robot_ (robot), comc_ (comc), joint_ (joint), reference_ (reference),
      mask_ (mask), nominalCase_ (false), jacobian_
      (3, robot->numberDof()-robot->extraConfigSpace().dimension())
    {
      if (mask[0] && mask[1] && mask[2])
        nominalCase_ = true;
      jacobian_.setZero ();
    }

    void RelativeCom::impl_compute (LiegroupElement& result,
				    ConfigurationIn_t argument)
      const throw ()
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      comc_->compute (Device::COM);
      const Transform3f& M = joint_->currentTransformation ();
      const vector3_t& x = comc_->com ();
      const matrix3_t& R = M.rotation ();
      const vector3_t& t = M.translation ();

      if (nominalCase_)
        result.vector () = R.transpose() * (x - t) - reference_;
      else {
        const vector3_t res ( R.transpose() * (x - t) - reference_);
        size_t index = 0;
        for (size_t i = 0; i < 3; ++i)
          if (mask_[i]) {
            result.vector () [index] = res [i];
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
      const matrix3_t& R (M.rotation ());
      const vector3_t& x (comc_->com ());
      const vector3_t& t (M.translation ());

      // Right part
      jacobian.rightCols (jacobian.cols () - Jjoint.cols ()).setZero ();
      // Left part
      // J = 0RTj ( Jcom + [ x - 0tj ]x 0Rj jJwj - 0Rj jJtj)
      jacobian_ = R.transpose() * Jcom;
      jacobian_.noalias() += (R.transpose() * R.colwise().cross(t-x)) * Jjoint.bottomRows<3>();

      if (nominalCase_) {
        jacobian.leftCols (Jjoint.cols ()).noalias() = jacobian_ - Jjoint.topRows<3>();
      } else {
        size_t index = 0;
        for (size_t i = 0; i < 3; ++i)
          if (mask_[i]) {
            jacobian.row(index).head(Jjoint.cols()) = jacobian_.row (i) - Jjoint.row(i);
            index++;
          }
      }
      hppDnum (info, "Jcom = " << std::endl << Jcom);
      hppDnum (info, "Jw = " << std::endl << Jjoint.bottomRows<3>());
      hppDnum (info, "Jv = " << std::endl << Jjoint.topRows<3>());
    }
  } // namespace constraints
} // namespace hpp
