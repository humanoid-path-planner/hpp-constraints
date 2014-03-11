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

#include <hpp/util/debug.hh>
#include <hpp/model/device.hh>
#include <hpp/model/fcl-to-eigen.hh>
#include <hpp/model/joint.hh>
#include <hpp/constraints/relative-orientation.hh>

namespace hpp {
  namespace constraints {

    extern void computeJlog (const double& theta, const vector_t& r,
			     eigen::matrix3_t& Jlog);

    RelativeOrientationPtr_t RelativeOrientation::create
    (const DevicePtr_t& robot, const JointPtr_t& joint1,
     const JointPtr_t& joint2, const matrix3_t& reference, bool ignoreZ)
    {
      RelativeOrientation* ptr =
	new RelativeOrientation (robot, joint1, joint2, reference, ignoreZ);
      RelativeOrientationShPtr shPtr (ptr);
      return shPtr;
    }

    RelativeOrientation::RelativeOrientation
    (const DevicePtr_t& robot, const JointPtr_t& joint1,
     const JointPtr_t& joint2, const matrix3_t& reference, bool ignoreZ) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
		ignoreZ ? 2 : 3, "RelativeOrientation"),
      robot_ (robot), joint1_ (joint1), joint2_ (joint2),
      reference_ (reference), ignoreZ_ (ignoreZ), r_ (3),
      Jlog_ (), jacobian_ (3, robot->numberDof ())
    {
    }

    void RelativeOrientation::computeError (vector_t& result,
					    const vector_t& argument,
					    double& theta, bool ignoreZ) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      const Transform3f& M1 = joint1_->currentTransformation ();
      const Transform3f& M2 = joint2_->currentTransformation ();
      const fcl::Matrix3f& R1 (M1.getRotation ());
      fcl::Matrix3f R2T (M2.getRotation ()); R2T.transpose ();
      Rerror_ = R2T * R1 * reference_;
      double tr = Rerror_ (0, 0) + Rerror_ (1, 1) + Rerror_ (2, 2);
      if (tr > 3) tr = 3;
      if (tr < -1) tr = -1;
      theta = acos ((tr - 1)/2);
      assert (theta == theta);
      if (theta > 1e-6) {
	result [0] = theta*(Rerror_ (2, 1) - Rerror_ (1, 2))/(2*sin(theta));
	result [1] = theta*(Rerror_ (0, 2) - Rerror_ (2, 0))/(2*sin(theta));
	if (!ignoreZ) {
	  result [2] = theta*(Rerror_ (1, 0) - Rerror_ (0, 1))/(2*sin(theta));
	}
      }
      else {
	result [0] = (Rerror_ (2, 1) - Rerror_ (1, 2))/2;
	result [1] = (Rerror_ (0, 2) - Rerror_ (2, 0))/2;
	if (!ignoreZ) {
	  result [2] = (Rerror_ (1, 0) - Rerror_ (0, 1))/2;
	}
      }
    }

    void RelativeOrientation::impl_compute (vector_t& result,
					    const vector_t& argument)
      const throw ()
    {
      double theta;
      computeError (result, argument, theta, ignoreZ_);
    }

    void RelativeOrientation::impl_jacobian (matrix_t &jacobian,
					     const vector_t &arg)
      const throw ()
    {
      const Transform3f& M2 = joint2_->currentTransformation ();
      fcl::Matrix3f R2T (M2.getRotation ()); R2T.transpose ();
      // Compute vector r
      double theta;
      computeError (r_, arg, theta, false);
      assert (theta >= 0);
      if (theta < 1e-3) {
	Jlog_.setIdentity ();
      } else {
	computeJlog (theta, r_, Jlog_);
      }
      const JointJacobian_t& Jjoint1 (joint1_->jacobian ());
      const JointJacobian_t& Jjoint2 (joint2_->jacobian ());
      hppDout (info, "Jlog_ = " << std::endl << Jlog_);
      hppDout (info, "Jjoint1 = " << std::endl << Jjoint1);
      hppDout (info, "Jjoint2 = " << std::endl << Jjoint2);
      if (ignoreZ_) {
	jacobian_ = Jlog_ * R2T *
	  (Jjoint1.bottomRows (3) - Jjoint2.bottomRows (3));
	jacobian = jacobian_.topRows (2);
      }
      else {
	hppDout (info, "Jlog_ * R2T = " << Jlog_ * R2T);
	jacobian = Jlog_ * R2T *
	  (Jjoint1.bottomRows (3) - Jjoint2.bottomRows (3));
      }
    }
  } // namespace constraints
} // namespace hpp

