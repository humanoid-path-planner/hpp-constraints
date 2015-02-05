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
        const JointPtr_t& jointRef, std::vector <bool> mask)
    {
      CenterOfMassComputationPtr_t comc =
        CenterOfMassComputation::create (robot);
      comc->add (robot->rootJoint ());
      comc->computeMass ();
      return create (name, robot, comc, jointL, jointR,
          pointL, pointR, jointRef, mask);
    }

    ComBetweenFeetPtr_t ComBetweenFeet::create (
        const std::string& name, const DevicePtr_t& robot,
        const CenterOfMassComputationPtr_t& comc,
        const JointPtr_t& jointL, const JointPtr_t& jointR,
        const vector3_t   pointL, const vector3_t   pointR,
        const JointPtr_t& jointRef, std::vector <bool> mask)
    {
      ComBetweenFeet* ptr = new ComBetweenFeet
        (name, robot, comc, jointL, jointR, pointL, pointR, jointRef, mask);
      ComBetweenFeetPtr_t shPtr (ptr);
      return shPtr;
    }

    ComBetweenFeet::ComBetweenFeet (
        const std::string& name, const DevicePtr_t& robot,
        const CenterOfMassComputationPtr_t& comc,
        const JointPtr_t& jointL, const JointPtr_t& jointR,
        const vector3_t   pointL, const vector3_t   pointR,
        const JointPtr_t& jointRef,
        std::vector <bool> mask) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
          size (mask), name),
      robot_ (robot), com_ (comc), left_ (jointL, pointL),
      right_ (jointR, pointR), jointRef_ (jointRef),
      mask_ (mask),
      result_ (3), jacobian_ (3, robot->numberDof ())
    {
      cross_.setZero ();
      jacobian_.setZero ();
      u_ = right_ - left_;
      xmxl_ = com_ - left_;
      xmxr_ = com_ - right_;
      ecrossu_ = (com_ - ((left_ + right_) * 0.5))^(u_);
      expr_ = RotationMultiply <ECrossU_t> (jointRef_, ecrossu_, true);
    }

    void ComBetweenFeet::impl_compute (vectorOut_t result,
        ConfigurationIn_t argument)
      const throw ()
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      size_t index = 0;
      u_.computeValue ();
      if (mask_[0]) {
        expr_.computeValue ();
        result[index++] = expr_.value ()[0];
      }
      if (mask_[1]) {
        xmxr_.computeValue ();
        result[index++] = - xmxr_.value().dot(u_.value ());
      }
      if (mask_[2]) {
        xmxl_.computeValue ();
        result[index  ] =   xmxl_.value().dot(u_.value ());
      }

      // com_->compute (Device::COM);
      // left_ .computeGlobal ();
      // right_.computeGlobal ();
      // const Transform3f& Mref = jointRef_->currentTransformation ();
      // eigen::matrix3_t RT; convert (Mref.getRotation (), RT);
      // RT.transposeInPlace ();
      // eigen::vector3_t x; convert (comc_->com (), x);
      // const eigen::vector3_t&
        // xl =  left_.global (), xr = right_.global (),
        // e = x  - (xl + xr)*0.5, u = xr - xl;
      // computeCrossMatrix (e, cross_);
      // size_t index = 0;
      // if (mask_[0]) result[index++] = (RT.row (0)*cross_*u);
      // if (mask_[1]) result[index++] = - (x - xr).dot(u);
      // if (mask_[2]) result[index  ] =   (x - xl).dot(u);
    }

    void ComBetweenFeet::impl_jacobian (matrixOut_t jacobian,
        ConfigurationIn_t arg) const throw ()
    {
      robot_->currentConfiguration (arg);
      robot_->computeForwardKinematics ();
      size_t index = 0;
      u_.computeJacobian ();
      if (mask_[0]) {
        expr_.computeJacobian ();
        jacobian.row (index++).leftCols (jointRef_->jacobian ().cols ())
          = expr_.jacobian ().row (0);
      }
      if (mask_[1]) {
        xmxr_.computeJacobian ();
        jacobian.row (index++).leftCols (jointRef_->jacobian ().cols ())
          = - u_.value ().transpose () * xmxr_.jacobian ()
            - xmxr_.value ().transpose () * u_.jacobian ();
      }
      if (mask_[2]) {
        xmxl_.computeJacobian ();
        jacobian.row (index  ).leftCols (jointRef_->jacobian ().cols ())
          = - u_.value ().transpose () * xmxl_.jacobian ()
            - xmxl_.value ().transpose () * u_.jacobian ();
      }
      // comc_->compute (Device::JACOBIAN);
      // CrossProduct cp = left_ ^ rigth_;
      // left_ .computeGlobal (); left_ .computeJacobian ();
      // right_.computeGlobal (); right_.computeJacobian ();
      // const ComJacobian_t& Jcom = comc_->jacobian ();
      // const JointJacobian_t& Jref (jointRef_->jacobian ());
      // eigen::vector3_t x; convert (comc_->com (), x);
      // const Transform3f& Mref = jointRef_->currentTransformation ();
      // eigen::matrix3_t RT; convert (Mref.getRotation (), RT);
      // RT.transposeInPlace ();

      // const eigen::vector3_t&
        // xl =  left_.global (), xr = right_.global (),
        // e = x  - (xl + xr) / 2, u = xr - xl;
      // eigen::matrix3_t ucross, eucross; ucross.setZero (); eucross.setZero ();
      // cross (u, ucross); cross (- ucross * e, eucross);
      // eigen::matrix3_t xcross; x.setZero ();
      // cross (c, xcross);
      // eigen::matrix3_t xmxrcross, xmxlcross; xmxrcross.setZero (); xmxlcross.setZero ();
      // cross (x - xr, xmxrcross); cross (x - xl, xmxlcross);
      // eigen::matrix3_t Rrprcross, Rlplcross; Rrprcross.setZero (); Rlplcross.setZero ();
      // cross (Mr.getRotation () * pointR_, Rrprcross);
      // cross (Ml.getRotation () * pointL_, Rlplcross);

      // size_t index = 0;
      // if (mask_[0])
        // jacobian.row (index++).leftCols (Jref.cols ()) = eigenRT.row(0) * (
            // eucross * Jref.bottomRows (3)
            // - ucross * Jcom
            // + xmxlcross * ( - Rrprcross * Jr.bottomRows (3) + Jr.topRows (3) )
            // - xmxrcross * ( - Rlplcross * Jl.bottomRows (3) + Jl.topRows (3) )
            // );
      // if (mask_[1])
      // jacobian.rightCols (jacobian.cols () - Jref.cols ()).setZero ();
    }
  } // namespace constraints
} // namespace hpp
