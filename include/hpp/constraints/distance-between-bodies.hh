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

#ifndef HPP_CONSTRAINTS_DISTANCE_BETWEEN_BODIES_HH
# define HPP_CONSTRAINTS_DISTANCE_BETWEEN_BODIES_HH

# include <pinocchio/multibody/geometry.hpp>

# include <hpp/pinocchio/collision-object.hh>
# include <hpp/pinocchio/liegroup-element.hh>

# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/differentiable-function.hh>

namespace hpp {
  namespace constraints {
    /// Distance between two sets of objects
    ///
    /// This function maps to a configuration of a robot, the distance
    ///   \li either between objects of a joints and objects of another joint,
    ///   \li or objects of a joint with a list of fixed objects.
    ///
    /// The above type of distance is determined by the method "create" called.
    class HPP_CONSTRAINTS_DLLAPI DistanceBetweenBodies :
      public DifferentiableFunction
    {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      
      /// Create instance and return shared pointer
      ///
      /// \param name name of the constraint,
      /// \param robot robot that own the bodies,
      /// \param joint1 joint that holds the first body,
      /// \param joint2 joint that holds the second body.
      static DistanceBetweenBodiesPtr_t create (const std::string& name,
						const DevicePtr_t& robot,
						const JointPtr_t& joint1,
						const JointPtr_t& joint2);
      
      /// Create instance and return shared pointer
      ///
      /// \param name name of the constraint,
      /// \param robot robot that own the bodies,
      /// \param joint joint that holds the body,
      /// \param objects list of fixed objects in the environment.
      static DistanceBetweenBodiesPtr_t create (const std::string& name,
						const DevicePtr_t& robot,
						const JointPtr_t& joint,
						const std::vector<CollisionObjectPtr_t>& objects);

      virtual ~DistanceBetweenBodies () throw () {}

    protected:
      /// Protected constructor
      ///
      /// \param name name of the constraint,
      /// \param robot robot that own the bodies,
      /// \param joint1 joint that holds the first body,
      /// \param joint2 joint that holds the second body.
      DistanceBetweenBodies (const std::string& name, const DevicePtr_t& robot,
			     const JointPtr_t& joint1,
			     const JointPtr_t& joint2);

      /// Protected constructor
      ///
      /// \param name name of the constraint,
      /// \param robot robot that own the bodies,
      /// \param joint joint that holds the body,
      /// \param objects list of fixed objects in the environment.
      DistanceBetweenBodies (const std::string& name, const DevicePtr_t& robot,
			     const JointPtr_t& joint,
			     const std::vector<CollisionObjectPtr_t>& objects);

      virtual void impl_compute (LiegroupElementRef result,
				 ConfigurationIn_t argument) const throw ();
      virtual void impl_jacobian (matrixOut_t jacobian,
				  ConfigurationIn_t arg) const throw ();
    private:
      typedef ::pinocchio::GeometryData GeometryData;

      DevicePtr_t robot_;
      JointPtr_t joint1_;
      JointPtr_t joint2_;
      mutable GeometryData data_;
      mutable std::size_t minIndex_;
      mutable Configuration_t latestArgument_;
      mutable LiegroupElement latestResult_;
    }; // class DistanceBetweenBodies
  } // namespace constraints
} // namespace hpp

#endif //HPP_CONSTRAINTS_DISTANCE_BETWEEN_BODIES_HH
