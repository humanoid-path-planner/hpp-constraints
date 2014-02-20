//
// Copyright (c) 2014 CNRS
// Authors: Florent Lamiraux
//
//
// This file is part of hpp-wholebody-step.
// hpp-wholebody-step-planner is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-wholebody-step-planner is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-wholebody-step-planner. If not, see
// <http://www.gnu.org/licenses/>.

#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include <hpp/constraints/orientation.hh>

namespace hpp {
  namespace constraints {

    static vector3d zero3d (0, 0, 0);

    void extractRotation (const matrix4d& M, matrix3d& R)
    {
      R (0,0) = M(0,0); R (0,1) = M(0,1); R (0,2) = M(0,2);
      R (1,0) = M(1,0); R (1,1) = M(1,1); R (1,2) = M(1,2);
      R (2,0) = M(2,0); R (2,1) = M(2,1); R (2,2) = M(2,2);
    }

    void computeJlog (const double& theta, const vector3d& r, matrix3d& Jlog)
    {
      if (theta < 1e-6)
	Jlog.setIdentity ();
      else {
	Jlog.setZero ();
	// Jlog = alpha I
	double alpha = .5*theta*sin(theta)/(1-cos (theta));
	Jlog (0,0) = Jlog (1,1) = Jlog (2,2) = alpha;
	// Jlog += -r_{\times}/2
	Jlog (0,1) = .5*r [2]; Jlog (1,0) = -.5*r [2];
	Jlog (0,2) = -.5*r [1]; Jlog (2,0) = .5*r [1];
	Jlog (1,2) = .5*r [0]; Jlog (2,1) = -.5*r [0];
	alpha = 1/(theta*theta) - sin(theta)/(2*theta*(1-cos(theta)));
	Jlog += alpha * r * r.transpose ();
      }
    }

    OrientationPtr_t Orientation::create (const DevicePtr_t& robot,
					  const JointPtr_t& joint,
					  const matrix3d& reference,
					  bool ignoreZ)
    {
      Orientation* ptr = new Orientation (robot, joint, reference, ignoreZ);
      OrientationShPtr shPtr (ptr);
      return shPtr;
    }
    Orientation::Orientation (const DevicePtr_t& robot, const JointPtr_t& joint,
			      const matrix3d& reference, bool ignoreZ) :
      parent_t (robot->numberDof (), ignoreZ ? 2 : 3, "Orientation"),
      robot_ (robot), joint_ (joint), reference_ (reference),
      ignoreZ_ (ignoreZ), r_ (3), Jlog_ (),
      jacobian_ (3, robot->numberDof ()), Jjoint_ (6, robot->numberDof ())
    {
    }

    void Orientation::computeError (result_t& result,
				    const argument_t& argument,
				    double& theta, bool ignoreZ) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      const matrix4d& M = joint_->currentTransformation ();
      extractRotation (M, R_);
      Rerror_ = R_.transpose () * reference_;
      double tr = Rerror_ (0, 0) + Rerror_ (1, 1) + Rerror_ (2, 2);
      if (tr > 1) tr = 1;
      if (tr < -1) tr = -1;
      theta = acos ((tr - 1)/2);
      result [0] = theta*(Rerror_ (2, 1) - Rerror_ (1, 2))/(2*sin(theta));
      result [1] = theta*(Rerror_ (0, 2) - Rerror_ (2, 0))/(2*sin(theta));
      if (!ignoreZ) {
	result [2] = theta*(Rerror_ (1, 0) - Rerror_ (0, 1))/(2*sin(theta));
      }
    }

    void Orientation::impl_compute (result_t& result,
				    const argument_t& argument) const throw ()
    {
      double theta;
      computeError (result, argument, theta, ignoreZ_);
    }
    void Orientation::impl_jacobian (jacobian_t &jacobian,
				     const argument_t &arg) const throw ()
    {
      // Compute vector r
      double theta;
      computeError (r_, arg, theta, false);
      assert (theta >= 0);
      if (theta < 1e-3) {
	Jlog_.setIdentity ();
      } else {
	computeJlog (theta, r_, Jlog_);
	robot_->getJacobian (*(robot_->rootJoint ()), *(joint_->dynamic ()),
			     zero3d, Jjoint_);
	jacobian_ = -Jlog_ * R_.transpose () * Jjoint_.bottomRows (3);
	if (ignoreZ_)
	  jacobian = jacobian_.topRows (2);
	else
	  jacobian = jacobian_;
      }
    }

    void Orientation::impl_gradient (gradient_t &gradient,
				     const argument_t &argument,
				     size_type functionId) const throw ()
    {
      matrix_t jacobian (outputSize (), inputSize ());
      impl_jacobian (jacobian, argument);
      gradient = jacobian.row (functionId);
    }

  } // namespace constraints
} // namespace hpp
