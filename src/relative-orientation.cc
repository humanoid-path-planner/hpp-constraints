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
#include <hpp/model/joint.hh>
#include <hpp/constraints/relative-orientation.hh>

namespace hpp {
  namespace constraints {

    static vector3d zero3d (0, 0, 0);
    extern void extractRotation (const matrix4d& M, matrix3d& R);
    extern void computeJlog (const double& theta, const vector3d& r,
			     matrix3d& Jlog);

    RelativeOrientationPtr_t RelativeOrientation::create
    (const DevicePtr_t& robot, const JointPtr_t& joint1,
     const JointPtr_t& joint2, const matrix3d& reference, bool ignoreZ)
    {
      RelativeOrientation* ptr =
	new RelativeOrientation (robot, joint1, joint2, reference, ignoreZ);
      RelativeOrientationShPtr shPtr (ptr);
      return shPtr;
    }

    RelativeOrientation::RelativeOrientation
    (const DevicePtr_t& robot, const JointPtr_t& joint1,
     const JointPtr_t& joint2, const matrix3d& reference, bool ignoreZ) :
      parent_t (robot->numberDof (), ignoreZ ? 2 : 3, "RelativeOrientation"),
      robot_ (robot), joint1_ (joint1), joint2_ (joint2),
      reference_ (reference), ignoreZ_ (ignoreZ), r_ (3),
      Jlog_ (), jacobian_ (3, robot->numberDof ()),
      Jjoint1_ (6, robot->numberDof ()), Jjoint2_ (6, robot->numberDof ())
    {
    }

    void RelativeOrientation::computeError (result_t& result,
					    const argument_t& argument,
					    double& theta, bool ignoreZ) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      const matrix4d& M1 = joint1_->currentTransformation ();
      const matrix4d& M2 = joint2_->currentTransformation ();
      extractRotation (M1, R1_);
      extractRotation (M2, R2_);
      Rerror_ = R2_.transpose () * R1_ * reference_;
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

    void RelativeOrientation::impl_compute (result_t& result,
					    const argument_t& argument)
      const throw ()
    {
      double theta;
      computeError (result, argument, theta, ignoreZ_);
    }

    void RelativeOrientation::impl_jacobian (jacobian_t &jacobian,
					     const argument_t &arg)
      const throw ()
    {
      // Compute vector r
      double theta;
      computeError (r_, arg, theta, false);
      assert (theta >= 0);
      if (theta < 1e-3) {
	Jlog_.setIdentity ();
      } else {
	computeJlog (theta, r_, Jlog_);
      }
      robot_->getJacobian (*(robot_->rootJoint ()), *(joint1_->dynamic ()),
			   zero3d, Jjoint1_);
      robot_->getJacobian (*(robot_->rootJoint ()), *(joint2_->dynamic ()),
			   zero3d, Jjoint2_);
      hppDout (info, "Jlog_ = " << std::endl << Jlog_);
      hppDout (info, "Jjoint1_ = " << std::endl << Jjoint1_);
      hppDout (info, "Jjoint2_ = " << std::endl << Jjoint2_);
      jacobian_ = Jlog_ * R2_.transpose () *
	(Jjoint1_.bottomRows (3) - Jjoint2_.bottomRows (3));
      if (ignoreZ_)
	jacobian = jacobian_.topRows (2);
      else
	jacobian = jacobian_;
    }

    void RelativeOrientation::impl_gradient (gradient_t &gradient,
					     const argument_t &argument,
					     size_type functionId)
      const throw ()
    {
      matrix_t jacobian (outputSize (), inputSize ());
      impl_jacobian (jacobian, argument);
      gradient = jacobian.row (functionId);
    }

  } // namespace constraints
} // namespace hpp

