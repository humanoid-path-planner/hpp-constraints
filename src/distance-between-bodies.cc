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
#include <hpp/pinocchio/joint-collection.hh>

namespace hpp {
  namespace constraints {
    typedef std::vector<CollisionObjectPtr_t> CollisionObjects_t;

    void initGeomData(const pinocchio::GeomModel& model, pinocchio::GeomData& data,
        const pinocchio::BodyPtr_t& body,
        const CollisionObjects_t& objects)
    {
      // Deactivate all collision pairs.
      for (std::size_t i = 0; i < model.collisionPairs.size(); ++i)
        data.activateCollisionPair(i, false);
      // Activate only the relevant ones.
      for (size_type i = 0; i < body->nbInnerObjects(); ++i) {
	CollisionObjectConstPtr_t obj1 (body->innerObjectAt(i));
	for (std::size_t j = 0; j < objects.size(); ++j) {
	  CollisionObjectConstPtr_t obj2 (objects[j]);
          std::size_t idx = model.findCollisionPair(
              ::pinocchio::CollisionPair (obj1->indexInModel(), obj2->indexInModel())
              );
          if (idx < model.collisionPairs.size())
            data.activateCollisionPair(idx);
          else
            throw std::invalid_argument("Collision pair not found");
	}
      }
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
     const JointPtr_t& joint, const CollisionObjects_t& objects)
    {
      DistanceBetweenBodies* ptr = new DistanceBetweenBodies
	(name, robot, joint, objects);
      DistanceBetweenBodiesPtr_t shPtr (ptr);
      return shPtr;
    }

    DistanceBetweenBodies::DistanceBetweenBodies
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const JointPtr_t& joint2) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
                              LiegroupSpace::R1 (), name), robot_ (robot),
      joint1_ (joint1), joint2_ (joint2), data_ (robot->geomModel()),
      latestResult_ (outputSpace ())
    {
      pinocchio::BodyPtr_t body2 (joint2_->linkedBody());
      CollisionObjects_t objects2 (body2->nbInnerObjects());
      for (std::size_t j = 0; j < objects2.size(); ++j)
        objects2[j] = body2->innerObjectAt (j);
      initGeomData(robot_->geomModel(), data_, joint1_->linkedBody(), objects2);
    }

    DistanceBetweenBodies::DistanceBetweenBodies
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint, const CollisionObjects_t& objects):
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
                              LiegroupSpace::R1 (), name),
      robot_ (robot), joint1_ (joint), joint2_ (), data_ (robot->geomModel()),
      latestResult_ (outputSpace ())
    {
      initGeomData(robot_->geomModel(), data_, joint1_->linkedBody(), objects);
    }

    void DistanceBetweenBodies::impl_compute
    (LiegroupElementRef result, ConfigurationIn_t argument) const throw ()
    {
      if ((argument.rows () == latestArgument_.rows ()) &&
	  (argument == latestArgument_)) {
	result = latestResult_;
	return;
      }
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      ::pinocchio::updateGeometryPlacements(robot_->model(), robot_->data(), robot_->geomModel(), data_);
      minIndex_ = ::pinocchio::computeDistances(robot_->geomModel(), data_);
      result.vector () [0] = data_.distanceResults[minIndex_].min_distance;
      latestArgument_ = argument;
      latestResult_ = result;
    }

    void DistanceBetweenBodies::impl_jacobian
    (matrixOut_t jacobian, ConfigurationIn_t arg) const throw ()
    {
      LiegroupElement dist (outputSpace ());
      impl_compute (dist, arg);
      const JointJacobian_t& J1 (joint1_->jacobian());
      const Transform3f& M1 (joint1_->currentTransformation());
      const matrix3_t& R1 (M1.rotation());
      vector3_t point1 (data_.distanceResults[minIndex_].nearest_points[0]);
      vector3_t point2 (data_.distanceResults[minIndex_].nearest_points[1]);
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
                  ) / dist.vector () [0];
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
	jacobian.noalias() -= tmp2/dist.vector () [0];
      }
    }
  } // namespace constraints
} // namespace hpp
