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

#ifndef HPP_CONSTRAINTS_DISTANCE_BETWEEN_POINTS_IN_BODIES_HH
# define HPP_CONSTRAINTS_DISTANCE_BETWEEN_POINTS_IN_BODIES_HH

# include <hpp/constraints/differentiable-function.hh>

namespace hpp {
  namespace constraints {
    /// Distance between two sets of objects
    ///
    /// This function maps to a configuration of a robot, the distance
    ///   \li either between two points in two joints
    ///   \li or between a point in a joint and a point in the environment
    ///
    /// The type of distance above is determined by the method "create" called.
    class HPP_CONSTRAINTS_DLLAPI DistanceBetweenPointsInBodies :
      public DifferentiableFunction
    {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW

      /// Create instance and return shared pointer
      ///
      /// \param name name of the constraint,
      /// \param robot robot that own the bodies,
      /// \param joint1 joint that holds the first point,
      /// \param joint2 joint that holds the second point,
      /// \param point1 point in frame of joint 1,
      /// \param point2 point in frame of joint 2.
      static DistanceBetweenPointsInBodiesPtr_t create
	(const std::string& name, const DevicePtr_t& robot,
	 const JointPtr_t& joint1, const JointPtr_t& joint2,
	 const vector3_t& point1, const vector3_t& point2);

      /// Create instance and return shared pointer
      ///
      /// \param name name of the constraint,
      /// \param robot robot that own the bodies,
      /// \param joint1 joint that holds the first point,
      /// \param point1 point in frame of joint 1,
      /// \param point2 point in frame of joint 2.
      static DistanceBetweenPointsInBodiesPtr_t create
	(const std::string& name, const DevicePtr_t& robot,
	 const JointPtr_t& joint1, const vector3_t& point1,
	 const vector3_t& point2);

      virtual ~DistanceBetweenPointsInBodies () throw () {}

    protected:
      /// Protected constructor
      ///
      /// \param name name of the constraint,
      /// \param robot robot that own the bodies,
      /// \param joint1 joint that holds the first body,
      /// \param joint2 joint that holds the second body.
      DistanceBetweenPointsInBodies (const std::string& name,
				     const DevicePtr_t& robot,
				     const JointPtr_t& joint1,
				     const JointPtr_t& joint2,
				     const vector3_t& point1,
				     const vector3_t& point2);

      /// Protected constructor
      ///
      /// \param name name of the constraint,
      /// \param robot robot that own the bodies,
      /// \param joint1 joint that holds the first body,
      /// \param joint2 joint that holds the second body.
      DistanceBetweenPointsInBodies (const std::string& name,
				     const DevicePtr_t& robot,
				     const JointPtr_t& joint1,
				     const vector3_t& point1,
				     const vector3_t& point2);

      virtual void impl_compute (vectorOut_t result,
				 ConfigurationIn_t argument) const throw ();
      virtual void impl_jacobian (matrixOut_t jacobian,
				  ConfigurationIn_t arg) const throw ();
    private:
      DevicePtr_t robot_;
      JointPtr_t joint1_;
      JointPtr_t joint2_;
      vector3_t point1_;
      vector3_t point2_;
      mutable vector3_t global1_;
      mutable vector3_t global2_;
      mutable Configuration_t latestArgument_;
      mutable vector_t latestResult_;
    }; // class DistanceBetweenPointsInBodies
  } // namespace constraints
} // namespace hpp

#endif //HPP_CONSTRAINTS_DISTANCE_BETWEEN_POINTS_IN_BODIES_HH
