//
// Copyright (c) 2015 CNRS
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

#include <hpp/constraints/distance-between-points-in-bodies.hh>

#include <pinocchio/spatial/se3.hpp>

#include <hpp/pinocchio/collision-object.hh>
#include <hpp/pinocchio/body.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>

namespace hpp {
  namespace constraints {

    using hpp::pinocchio::LiegroupElement;
    using hpp::pinocchio::LiegroupSpace;

    static void cross (const vector3_t& v, eigen::matrix3_t& m)
    {
      m (0,1) = -v [2]; m (1,0) =  v [2];
      m (0,2) =  v [1]; m (2,0) = -v [1];
      m (1,2) = -v [0]; m (2,1) =  v [0];
      m (0,0) = m (1,1) = m (2,2) = 0;
    }

    DistanceBetweenPointsInBodiesPtr_t DistanceBetweenPointsInBodies::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const JointPtr_t& joint2,
     const vector3_t& point1, const vector3_t& point2)
    {
      DistanceBetweenPointsInBodies* ptr = new DistanceBetweenPointsInBodies
	(name, robot, joint1, joint2, point1, point2);
      DistanceBetweenPointsInBodiesPtr_t shPtr (ptr);
      return shPtr;
    }

    DistanceBetweenPointsInBodiesPtr_t DistanceBetweenPointsInBodies::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const vector3_t& point1, const vector3_t& point2)
    {
      DistanceBetweenPointsInBodies* ptr = new DistanceBetweenPointsInBodies
	(name, robot, joint1, point1, point2);
      DistanceBetweenPointsInBodiesPtr_t shPtr (ptr);
      return shPtr;
    }

    DistanceBetweenPointsInBodies::DistanceBetweenPointsInBodies
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const JointPtr_t& joint2,
     const vector3_t& point1, const vector3_t& point2) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
                              LiegroupSpace::R1 (), name), robot_ (robot),
      joint1_ (joint1), joint2_ (joint2), point1_ (point1), point2_ (point2),
      latestResult_ (outputSpace ())
    {
      assert (joint1);
      global2_ = point2;
    }

    DistanceBetweenPointsInBodies::DistanceBetweenPointsInBodies
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const vector3_t& point1, const vector3_t& point2)
      : DifferentiableFunction (robot->configSize (), robot->numberDof (),
                                LiegroupSpace::R1 (), name), robot_ (robot),
        joint1_ (joint1), joint2_ (), point1_ (point1), point2_ (point2),
        latestResult_ (outputSpace ())
    {
      assert (joint1);
      global2_ = point2;
    }

    void DistanceBetweenPointsInBodies::impl_compute
    (LiegroupElement& result, ConfigurationIn_t argument) const throw ()
    {
      if ((argument.rows () == latestArgument_.rows ()) &&
	  (argument == latestArgument_)) {
	result = latestResult_;
	return;
      }
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      global1_ = joint1_->currentTransformation ().act (point1_);
      if (joint2_) {
	global2_ = joint2_->currentTransformation ().act (point2_);
        result.vector () [0] = (global2_ - global1_).norm ();
      } else {
        result.vector () [0] = (           global1_).norm ();
      }

      latestArgument_ = argument;
      latestResult_ = result;
    }

    void DistanceBetweenPointsInBodies::impl_jacobian
    (matrixOut_t jacobian, ConfigurationIn_t arg) const throw ()
    {
      LiegroupElement dist (outputSpace ());
      impl_compute (dist, arg);
      const JointJacobian_t& J1 (joint1_->jacobian());
      const Transform3f& M1 (joint1_->currentTransformation());
      const matrix3_t& R1 (M1.rotation());

      // P1 - P2
      vector3_t P1_minus_P2 (global1_ - global2_);
      // P1 - t1
      vector3_t P1_minus_t1 (global1_ - M1.translation ());

      // FIXME Remove me
      eigen::matrix3_t P1_minus_t1_cross; cross(P1_minus_t1, P1_minus_t1_cross);
      assert (R1.colwise().cross(P1_minus_t1).isApprox(- P1_minus_t1_cross * R1));

      //        T (                              )
      // (P1-P2)  ( J    -   [P1 - t1]  J        )
      //          (  1 [0:3]          x  1 [3:6] )
      matrix_t tmp1
	(P1_minus_P2.transpose () * R1 * J1.topRows (3)
         + P1_minus_P2.transpose () * R1.colwise().cross(P1_minus_t1) * J1.bottomRows (3));
      if (joint2_) {
        const JointJacobian_t& J2 (joint2_->jacobian());
        const Transform3f& M2 (joint2_->currentTransformation());
        const matrix3_t& R2 (M2.rotation());
	// P2 - t2
	vector3_t P2_minus_t2 (global2_ - M2.translation ());
	//        T (                              )
	// (P1-P2)  ( J    -   [P1 - t1]  J        )
	//          (  2 [0:3]          x  2 [3:6] )
	matrix_t tmp2
	  (P1_minus_P2.transpose () * R2 * J2.topRows (3)
           + P1_minus_P2.transpose () * R2.colwise().cross(P2_minus_t2) * J2.bottomRows (3));
	jacobian = (tmp1 - tmp2)/dist.vector () [0];
      } else {
	jacobian = tmp1/dist.vector () [0];
      }
    }

  } // namespace constraints
} // namespace hpp
