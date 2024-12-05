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

#include <hpp/constraints/distance-between-points-in-bodies.hh>
#include <hpp/pinocchio/body.hh>
#include <hpp/pinocchio/collision-object.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <pinocchio/spatial/se3.hpp>

namespace hpp {
namespace constraints {

static void cross(const vector3_t& v, eigen::matrix3_t& m) {
  m(0, 1) = -v[2];
  m(1, 0) = v[2];
  m(0, 2) = v[1];
  m(2, 0) = -v[1];
  m(1, 2) = -v[0];
  m(2, 1) = v[0];
  m(0, 0) = m(1, 1) = m(2, 2) = 0;
}

DistanceBetweenPointsInBodiesPtr_t DistanceBetweenPointsInBodies::create(
    const std::string& name, const DevicePtr_t& robot, const JointPtr_t& joint1,
    const JointPtr_t& joint2, const vector3_t& point1,
    const vector3_t& point2) {
  DistanceBetweenPointsInBodies* ptr = new DistanceBetweenPointsInBodies(
      name, robot, joint1, joint2, point1, point2);
  DistanceBetweenPointsInBodiesPtr_t shPtr(ptr);
  return shPtr;
}

DistanceBetweenPointsInBodiesPtr_t DistanceBetweenPointsInBodies::create(
    const std::string& name, const DevicePtr_t& robot, const JointPtr_t& joint1,
    const vector3_t& point1, const vector3_t& point2) {
  DistanceBetweenPointsInBodies* ptr =
      new DistanceBetweenPointsInBodies(name, robot, joint1, point1, point2);
  DistanceBetweenPointsInBodiesPtr_t shPtr(ptr);
  return shPtr;
}

DistanceBetweenPointsInBodies::DistanceBetweenPointsInBodies(
    const std::string& name, const DevicePtr_t& robot, const JointPtr_t& joint1,
    const JointPtr_t& joint2, const vector3_t& point1, const vector3_t& point2)
    : DifferentiableFunction(robot->configSize(), robot->numberDof(),
                             LiegroupSpace::R1(), name),
      robot_(robot),
      joint1_(joint1),
      joint2_(joint2),
      point1_(point1),
      point2_(point2),
      latestResult_(outputSpace()) {
  assert(joint1);
  global2_ = point2;
}

DistanceBetweenPointsInBodies::DistanceBetweenPointsInBodies(
    const std::string& name, const DevicePtr_t& robot, const JointPtr_t& joint1,
    const vector3_t& point1, const vector3_t& point2)
    : DifferentiableFunction(robot->configSize(), robot->numberDof(),
                             LiegroupSpace::R1(), name),
      robot_(robot),
      joint1_(joint1),
      joint2_(),
      point1_(point1),
      point2_(point2),
      latestResult_(outputSpace()) {
  assert(joint1);
  global2_ = point2;
}

void DistanceBetweenPointsInBodies::impl_compute(
    LiegroupElementRef result, ConfigurationIn_t argument) const {
  if ((argument.rows() == latestArgument_.rows()) &&
      (argument == latestArgument_)) {
    result = latestResult_;
    return;
  }
  robot_->currentConfiguration(argument);
  robot_->computeForwardKinematics(pinocchio::JOINT_POSITION |
                                   pinocchio::JACOBIAN);
  global1_ = joint1_->currentTransformation().act(point1_);
  if (joint2_) {
    global2_ = joint2_->currentTransformation().act(point2_);
    result.vector()[0] = (global2_ - global1_).norm();
  } else {
    result.vector()[0] = (global1_).norm();
  }

  latestArgument_ = argument;
  latestResult_ = result;
}

void DistanceBetweenPointsInBodies::impl_jacobian(matrixOut_t jacobian,
                                                  ConfigurationIn_t arg) const {
  LiegroupElement dist(outputSpace());
  impl_compute(dist, arg);
  const JointJacobian_t& J1(joint1_->jacobian());
  const Transform3s& M1(joint1_->currentTransformation());
  const matrix3_t& R1(M1.rotation());

  // P1 - P2
  vector3_t P1_minus_P2(global1_ - global2_);
  // P1 - t1
  vector3_t P1_minus_t1(global1_ - M1.translation());

  // FIXME Remove me
  eigen::matrix3_t P1_minus_t1_cross;
  cross(P1_minus_t1, P1_minus_t1_cross);
  assert(R1.colwise().cross(P1_minus_t1).isApprox(-P1_minus_t1_cross * R1));

  //        T (                              )
  // (P1-P2)  ( J    -   [P1 - t1]  J        )
  //          (  1 [0:3]          x  1 [3:6] )
  matrix_t tmp1(P1_minus_P2.transpose() * R1 * J1.topRows(3) +
                P1_minus_P2.transpose() * R1.colwise().cross(P1_minus_t1) *
                    J1.bottomRows(3));
  if (joint2_) {
    const JointJacobian_t& J2(joint2_->jacobian());
    const Transform3s& M2(joint2_->currentTransformation());
    const matrix3_t& R2(M2.rotation());
    // P2 - t2
    vector3_t P2_minus_t2(global2_ - M2.translation());
    //        T (                              )
    // (P1-P2)  ( J    -   [P1 - t1]  J        )
    //          (  2 [0:3]          x  2 [3:6] )
    matrix_t tmp2(P1_minus_P2.transpose() * R2 * J2.topRows(3) +
                  P1_minus_P2.transpose() * R2.colwise().cross(P2_minus_t2) *
                      J2.bottomRows(3));
    jacobian = (tmp1 - tmp2) / dist.vector()[0];
  } else {
    jacobian = tmp1 / dist.vector()[0];
  }
}

}  // namespace constraints
}  // namespace hpp
