//
// Copyright (c) 2015 CNRS
// Authors: Florent Lamiraux
//
//

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#ifndef HPP_CONSTRAINTS_DISTANCE_BETWEEN_POINTS_IN_BODIES_HH
# define HPP_CONSTRAINTS_DISTANCE_BETWEEN_POINTS_IN_BODIES_HH

# include <hpp/pinocchio/liegroup-element.hh>
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

      virtual ~DistanceBetweenPointsInBodies () {}

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
      DistanceBetweenPointsInBodies (const std::string& name,
				     const DevicePtr_t& robot,
				     const JointPtr_t& joint1,
				     const vector3_t& point1,
				     const vector3_t& point2);

      virtual void impl_compute (LiegroupElementRef result,
				 ConfigurationIn_t argument) const;
      virtual void impl_jacobian (matrixOut_t jacobian,
				  ConfigurationIn_t arg) const;

      bool isEqual(const DifferentiableFunction& other) const {
        const DistanceBetweenPointsInBodies& castother = dynamic_cast<const DistanceBetweenPointsInBodies&>(other);
        if (!DifferentiableFunction::isEqual(other))
          return false;
        
        if (robot_ != castother.robot_)
          return false;
        if (joint1_ != castother.joint1_)
          return false;
        if (joint2_ != castother.joint2_)
          return false;
        if (point1_ != castother.point1_)
          return false;
        if (point2_ != castother.point2_)
          return false;
        
        return true;
      }
    private:
      DevicePtr_t robot_;
      JointPtr_t joint1_;
      JointPtr_t joint2_;
      vector3_t point1_;
      vector3_t point2_;
      mutable vector3_t global1_;
      mutable vector3_t global2_;
      mutable Configuration_t latestArgument_;
      mutable LiegroupElement latestResult_;
    }; // class DistanceBetweenPointsInBodies
  } // namespace constraints
} // namespace hpp

#endif //HPP_CONSTRAINTS_DISTANCE_BETWEEN_POINTS_IN_BODIES_HH
