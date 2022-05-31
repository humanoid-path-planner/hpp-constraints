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

#ifndef HPP_CONSTRAINTS_DISTANCE_BETWEEN_BODIES_HH
#define HPP_CONSTRAINTS_DISTANCE_BETWEEN_BODIES_HH

#include <hpp/constraints/differentiable-function.hh>
#include <hpp/constraints/fwd.hh>
#include <hpp/pinocchio/collision-object.hh>
#include <hpp/pinocchio/liegroup-element.hh>
#include <pinocchio/multibody/geometry.hpp>

namespace hpp {
namespace constraints {
/// Distance between two sets of objects
///
/// This function maps to a configuration of a robot, the distance
///   \li either between objects of a joints and objects of another joint,
///   \li or objects of a joint with a list of fixed objects.
///
/// The above type of distance is determined by the method "create" called.
class HPP_CONSTRAINTS_DLLAPI DistanceBetweenBodies
    : public DifferentiableFunction {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /// Create instance and return shared pointer
  ///
  /// \param name name of the constraint,
  /// \param robot robot that own the bodies,
  /// \param joint1 joint that holds the first body,
  /// \param joint2 joint that holds the second body.
  static DistanceBetweenBodiesPtr_t create(const std::string& name,
                                           const DevicePtr_t& robot,
                                           const JointPtr_t& joint1,
                                           const JointPtr_t& joint2);

  /// Create instance and return shared pointer
  ///
  /// \param name name of the constraint,
  /// \param robot robot that own the bodies,
  /// \param joint joint that holds the body,
  /// \param objects list of fixed objects in the environment.
  static DistanceBetweenBodiesPtr_t create(
      const std::string& name, const DevicePtr_t& robot,
      const JointPtr_t& joint,
      const std::vector<CollisionObjectPtr_t>& objects);

  virtual ~DistanceBetweenBodies() {}

 protected:
  /// Protected constructor
  ///
  /// \param name name of the constraint,
  /// \param robot robot that own the bodies,
  /// \param joint1 joint that holds the first body,
  /// \param joint2 joint that holds the second body.
  DistanceBetweenBodies(const std::string& name, const DevicePtr_t& robot,
                        const JointPtr_t& joint1, const JointPtr_t& joint2);

  /// Protected constructor
  ///
  /// \param name name of the constraint,
  /// \param robot robot that own the bodies,
  /// \param joint joint that holds the body,
  /// \param objects list of fixed objects in the environment.
  DistanceBetweenBodies(const std::string& name, const DevicePtr_t& robot,
                        const JointPtr_t& joint,
                        const std::vector<CollisionObjectPtr_t>& objects);

  virtual void impl_compute(LiegroupElementRef result,
                            ConfigurationIn_t argument) const;
  virtual void impl_jacobian(matrixOut_t jacobian, ConfigurationIn_t arg) const;

  bool isEqual(const DifferentiableFunction& other) const {
    const DistanceBetweenBodies& castother =
        dynamic_cast<const DistanceBetweenBodies&>(other);
    if (!DifferentiableFunction::isEqual(other)) return false;

    if (robot_ != castother.robot_) return false;
    if (joint1_ != castother.joint1_) return false;
    if (joint2_ != castother.joint2_) return false;
    if (data_ != castother.data_) return false;

    return true;
  }

 private:
  typedef ::pinocchio::GeometryData GeometryData;

  DevicePtr_t robot_;
  JointPtr_t joint1_;
  JointPtr_t joint2_;
  mutable GeometryData data_;
  mutable std::size_t minIndex_;
  mutable Configuration_t latestArgument_;
  mutable LiegroupElement latestResult_;
};  // class DistanceBetweenBodies
}  // namespace constraints
}  // namespace hpp

#endif  // HPP_CONSTRAINTS_DISTANCE_BETWEEN_BODIES_HH
