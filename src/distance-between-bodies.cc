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

#include <hpp/constraints/distance-between-bodies.hh>

#include <pinocchio/algorithm/geometry.hpp>

#include <hpp/pinocchio/body.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>

namespace hpp {
  namespace constraints {

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
      joint2_ (joint2), objs1_ (joint1_->linkedBody ()->innerObjects()),
      objs2_ (joint2_->linkedBody ()->innerObjects ()),
      data_ (robot->geomModel())
    {
      initGeomData();
    }

    DistanceBetweenBodies::DistanceBetweenBodies
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint, const ObjectVector_t& objects) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (), 1,
			      name), robot_ (robot), joint1_ (joint),
      joint2_ (), objs1_ (joint1_->linkedBody ()->innerObjects ()),
      objs2_ (objects),
      data_ (robot->geomModel())
    {
      initGeomData();
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
      se3::updateGeometryPlacements(robot_->model(), robot_->data(), robot_->geomModel(), data_);
      minIndex_ = se3::computeDistances(data_);
      result [0] = data_.distance_results[minIndex_].distance();
      latestArgument_ = argument;
      latestResult_ = result;
    }

    void DistanceBetweenBodies::impl_jacobian
    (matrixOut_t jacobian, ConfigurationIn_t arg) const throw ()
    {
      vector_t dist; dist.resize (1);
      impl_compute (dist, arg);
      const JointJacobian_t& J1 (joint1_->jacobian());
      const Transform3f& M1 (joint1_->currentTransformation());
      const matrix3_t& R1 (M1.rotation());
      vector3_t point1 (data_.distance_results[minIndex_].closestPointInner());
      vector3_t point2 (data_.distance_results[minIndex_].closestPointOuter());
      // P1 - P2
      vector3_t P1_minus_P2 (point1 - point2);
      // P1 - t1
      vector3_t P1_minus_t1 (point1 - M1.translation ());
      //        T (                              )
      // (P1-P2)  ( J    -   [P1 - t1]  J        )
      //          (  1 [0:3]          x  1 [3:6] )
      jacobian = (
          P1_minus_P2.transpose () * R1 * J1.topRows<3>()
          + P1_minus_P2.transpose () * R1.colwise().cross(P1_minus_t1) * J1.bottomRows<3>()
          ) / dist[0];
      if (joint2_) {
        const JointJacobian_t& J2 (joint2_->jacobian());
        const Transform3f& M2 (joint2_->currentTransformation());
        const matrix3_t& R2 (M2.rotation());
	// P2 - t2
	vector3_t P2_minus_t2 (point2 - M2.translation ());
	//        T (                              )
	// (P1-P2)  ( J    -   [P1 - t1]  J        )
	//          (  2 [0:3]          x  2 [3:6] )
	matrix_t tmp2
	  (  P1_minus_P2.transpose () * R2 * J2.topRows<3>()
           + P1_minus_P2.transpose () * R2.colwise().cross(P2_minus_t2) * J2.bottomRows<3>());
	jacobian.noalias() -= tmp2/dist [0];
      }
    }

    void DistanceBetweenBodies::initGeomData()
    {
      using se3::GeometryModel;
      const GeometryModel& model = robot_->geomModel();
      // Deactivate all collision pairs.
      for (GeometryModel::Index i = 0; i < model.collisionPairs.size(); ++i)
        data_.activateCollisionPair(i, false);
      // Activate only the relevant ones.
      for (ObjectVector_t::const_iterator it1 = objs1_.begin ();
	   it1 != objs1_.end (); ++it1) {
	CollisionObjectConstPtr_t obj1 (*it1);
	for (ObjectVector_t::const_iterator it2 = objs2_.begin ();
	     it2 != objs2_.end (); ++it2) {
	  CollisionObjectConstPtr_t obj2 (*it2);
          GeometryModel::Index idx = model.findCollisionPair(
              se3::CollisionPair (obj1->indexInModel(), obj2->indexInModel())
              );
          if (idx < model.collisionPairs.size())
            data_.activateCollisionPair(idx);
          else
            throw std::invalid_argument("Collision pair not found");
	}
      }
    }
  } // namespace constraints
} // namespace hpp
