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

#include <hpp/fcl/distance.h>
#include <hpp/model/collision-object.hh>
#include <hpp/model/body.hh>
#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include <hpp/constraints/distance-between-bodies.hh>

namespace hpp {
  namespace constraints {

    static void cross (const fcl::Vec3f& v, eigen::matrix3_t& m)
    {
      m (0,1) = -v [2]; m (1,0) =  v [2];
      m (0,2) =  v [1]; m (2,0) = -v [1];
      m (1,2) = -v [0]; m (2,1) =  v [0];
      m (0,0) = m (1,1) = m (2,2) = 0;
    }

    static void fclToEigen (const fcl::Vec3f& v, eigen::vector3_t& res)
    {
      res [0] = v[0]; res [1] = v[1]; res [2] = v[2];
    }

    DistanceBetweenBodiesPtr_t DistanceBetweenBodies::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const JointPtr_t& joint2)
    {
      DistanceBetweenBodies* ptr = new DistanceBetweenBodies
	(name, robot, joint1, joint2);
      DistanceBetweenBodiesPtr_t shPtr (ptr);
      return shPtr;
    }

    DistanceBetweenBodiesPtr_t DistanceBetweenBodies::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint, const ObjectVector_t& objects)
    {
      DistanceBetweenBodies* ptr = new DistanceBetweenBodies
	(name, robot, joint, objects);
      DistanceBetweenBodiesPtr_t shPtr (ptr);
      return shPtr;
    }

    DistanceBetweenBodies::DistanceBetweenBodies
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const JointPtr_t& joint2) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (), 1,
			      name), robot_ (robot), joint1_ (joint1),
      joint2_ (joint2), objs1_ (joint1_->linkedBody ()->innerObjects
				(hpp::model::DISTANCE)),
      objs2_ (joint2_->linkedBody ()->innerObjects (hpp::model::DISTANCE))
    {
    }

    DistanceBetweenBodies::DistanceBetweenBodies
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint, const ObjectVector_t& objects) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (), 1,
			      name), robot_ (robot), joint1_ (joint),
      joint2_ (0x0), objs1_ (joint1_->linkedBody ()->innerObjects
			     (hpp::model::DISTANCE)),
      objs2_ (objects)
    {
    }

    void DistanceBetweenBodies::impl_compute
    (vectorOut_t result, ConfigurationIn_t argument) const throw ()
    {
      if ((argument.rows () == latestArgument_.rows ()) &&
	  (argument == latestArgument_)) {
	result = latestResult_;
	return;
      }
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      fcl::DistanceRequest distanceRequest (true, 0, 0, fcl::GST_INDEP);
      result [0] = std::numeric_limits <value_type>::infinity ();
      for (ObjectVector_t::const_iterator it1 = objs1_.begin ();
	   it1 != objs1_.end (); ++it1) {
	CollisionObjectPtr_t obj1 (*it1);
	for (ObjectVector_t::const_iterator it2 = objs2_.begin ();
	     it2 != objs2_.end (); ++it2) {
	  CollisionObjectPtr_t obj2 (*it2);
	  fcl::DistanceResult distanceResult;
	  fcl::distance (obj1->fcl ().get (), obj2->fcl ().get (),
			 distanceRequest, distanceResult);
	  if (distanceResult.min_distance < result [0]) {
	    result [0] = distanceResult.min_distance;
	    point1_ = distanceResult.nearest_points [0];
	    point2_ = distanceResult.nearest_points [1];
	  }
	}
      }
      latestArgument_ = argument;
      latestResult_ = result;
    }

    void DistanceBetweenBodies::impl_jacobian
    (matrixOut_t jacobian, ConfigurationIn_t arg) const throw ()
    {
      vector_t dist; dist.resize (1);
      impl_compute (dist, arg);
      // P1 - P2
      eigen::vector3_t P1_minus_P2;
      fclToEigen (point1_ - point2_, P1_minus_P2);
      // P1 - t1
      vector3_t P1_minus_t1
	(point1_ - joint1_->currentTransformation ().getTranslation ());
      // [P1 - t1]
      //          x
      eigen::matrix3_t P1_minus_t1_cross;
      cross (P1_minus_t1, P1_minus_t1_cross);
      //        T (                              )
      // (P1-P2)  ( J    -   [P1 - t1]  J        )
      //          (  1 [0:3]          x  1 [3:6] )
      matrix_t tmp1
	(P1_minus_P2.transpose () * joint1_->jacobian ().topRows (3) -
	 P1_minus_P2.transpose () * P1_minus_t1_cross *
	 joint1_->jacobian ().bottomRows (3));
      if (joint2_) {
	// P2 - t2
	vector3_t P2_minus_t2
	  (point2_ - joint2_->currentTransformation ().getTranslation ());
	// [P2 - t2]
	//          x
	eigen::matrix3_t P2_minus_t2_cross;
	cross (P2_minus_t2, P2_minus_t2_cross);
	//        T (                              )
	// (P1-P2)  ( J    -   [P1 - t1]  J        )
	//          (  2 [0:3]          x  2 [3:6] )
	matrix_t tmp2
	  (P1_minus_P2.transpose () * joint2_->jacobian ().topRows (3) -
	   P1_minus_P2.transpose () * P2_minus_t2_cross *
	   joint2_->jacobian ().bottomRows (3));
	jacobian = (tmp1 - tmp2)/dist [0];
      } else {
	jacobian = tmp1/dist [0];
      }
    }

  } // namespace constraints
} // namespace hpp
