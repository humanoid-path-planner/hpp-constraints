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
#include <hpp/constraints/orientation.hh>
#include <hpp/constraints/tool.hh>

namespace hpp {
  namespace constraints {
    RelativeOrientationPtr_t RelativeOrientation::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const JointPtr_t& joint2,
     const matrix3_t& reference, std::vector <bool> mask)
    {
      RelativeOrientation* ptr =
	new RelativeOrientation (name, robot, joint1, joint2, reference, mask);
      RelativeOrientationPtr_t shPtr (ptr);
      return shPtr;
    }

    RelativeOrientationPtr_t RelativeOrientation::create
    (const DevicePtr_t& robot, const JointPtr_t& joint1,
     const JointPtr_t& joint2, const matrix3_t& reference,
     std::vector <bool> mask)
    {
      return create ("RelativeOrientation", robot, joint1, joint2, reference, mask);
    }

    RelativeOrientation::RelativeOrientation
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const JointPtr_t& joint2,
     const matrix3_t& reference, std::vector <bool> mask) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
			      Orientation::size (mask), name),
      robot_ (robot), joint1_ (joint1), joint2_ (joint2),
      reference_ (reference), mask_ (mask), r_ (3),
      Jlog_ (), jacobian_ (3, robot->numberDof ())
    {
      jacobian_.setZero ();
    }

    void RelativeOrientation::computeError
    (vectorOut_t result, ConfigurationIn_t argument, double& theta,
     std::vector <bool> mask) const
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
      size_type index = 0;
      if (theta > 1e-6) {
	if (mask [0]) {
	  result [index] = theta*(Rerror_ (2, 1) -
				  Rerror_ (1, 2))/(2*sin(theta)); ++index;
	}
	if (mask [1]) {
	  result [index] = theta*(Rerror_ (0, 2) -
				  Rerror_ (2, 0))/(2*sin(theta)); ++index;
	}
	if (mask [2]) {
	  result [index] = theta*(Rerror_ (1, 0) -
				  Rerror_ (0, 1))/(2*sin(theta)); ++index;
	}
      }
      else {
	if (mask [0]) {
	  result [index] = (Rerror_ (2, 1) - Rerror_ (1, 2))/2; ++index;
	}
	if (mask [1]) {
	  result [index] = (Rerror_ (0, 2) - Rerror_ (2, 0))/2; ++index;
	}
	if (mask [2]) {
	  result [index] = (Rerror_ (1, 0) - Rerror_ (0, 1))/2; ++index;
	}
      }
    }

    void RelativeOrientation::impl_compute (vectorOut_t result,
					    ConfigurationIn_t argument)
      const throw ()
    {
      double theta;
      computeError (result, argument, theta, mask_);
    }

    void RelativeOrientation::impl_jacobian (matrixOut_t jacobian,
					     ConfigurationIn_t arg)
      const throw ()
    {
      const Transform3f& M2 = joint2_->currentTransformation ();
      fcl::Matrix3f R2T (M2.getRotation ()); R2T.transpose ();
      // Compute vector r
      double theta;
      computeError (r_, arg, theta, boost::assign::list_of (true)(true)(true));
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
      jacobian_.leftCols (Jjoint1.cols ()) = Jlog_ * R2T *
	(Jjoint1.bottomRows (3) - Jjoint2.bottomRows (3));
      size_type index = 0;
      if (mask_ [0]) {
	jacobian.row (index) = jacobian_.row (0); ++index;
      }
      if (mask_ [1]) {
	jacobian.row (index) = jacobian_.row (1); ++index;
      }
      if (mask_ [2]) {
	jacobian.row (index) = jacobian_.row (2);
      }
    }
  } // namespace constraints
} // namespace hpp

