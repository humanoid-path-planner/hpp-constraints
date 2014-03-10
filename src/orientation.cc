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

#include <hpp/model/device.hh>
#include <hpp/model/fcl-to-eigen.hh>
#include <hpp/model/joint.hh>
#include <hpp/constraints/orientation.hh>

namespace hpp {
  namespace constraints {

    static vector3_t zero3d (0, 0, 0);

    void computeJlog (const double& theta, const vector_t& r,
		      eigen::matrix3_t& Jlog)
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
					  const matrix3_t& reference,
					  bool ignoreZ)
    {
      Orientation* ptr = new Orientation (robot, joint, reference, ignoreZ);
      OrientationShPtr shPtr (ptr);
      return shPtr;
    }
    Orientation::Orientation (const DevicePtr_t& robot, const JointPtr_t& joint,
			      const matrix3_t& reference, bool ignoreZ) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
		ignoreZ ? 2 : 3, "Orientation"),
      robot_ (robot), joint_ (joint), reference_ (reference),
      ignoreZ_ (ignoreZ), r_ (3), Jlog_ (),
      jacobian_ (3, robot->numberDof ())
    {
    }

    void Orientation::computeError (vector_t& result,
				    const vector_t& argument,
				    double& theta, bool ignoreZ) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      const Transform3f& M = joint_->currentTransformation ();
      fcl::Matrix3f RT (M.getRotation ()); RT.transpose ();
      Rerror_ = RT * reference_;
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

    void Orientation::impl_compute (vector_t& result,
				    const vector_t& argument) const throw ()
    {
      double theta;
      computeError (result, argument, theta, ignoreZ_);
    }
    void Orientation::impl_jacobian (matrix_t &jacobian,
				     const vector_t &arg) const throw ()
    {
      const Transform3f& M = joint_->currentTransformation ();
      fcl::Matrix3f RT (M.getRotation ()); RT.transpose ();
      // Compute vector r
      double theta;
      computeError (r_, arg, theta, false);
      assert (theta >= 0);
      if (theta < 1e-3) {
	Jlog_.setIdentity ();
      } else {
	computeJlog (theta, r_, Jlog_);
	const JointJacobian_t& Jjoint (joint_->jacobian ());
	if (ignoreZ_) {
	  jacobian_ = -Jlog_ * RT * Jjoint.bottomRows (3);
	  jacobian = jacobian_.topRows (2);
	}
	else {
	  jacobian = -Jlog_ * RT * Jjoint.bottomRows (3);
	}
      }
    }

  } // namespace constraints
} // namespace hpp
