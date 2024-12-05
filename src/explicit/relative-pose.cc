// Copyright (c) 2018, LAAS-CNRS
// Authors: Florent Lamiraux
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

#include <boost/serialization/weak_ptr.hpp>
#include <hpp/constraints/explicit/relative-pose.hh>
#include <hpp/constraints/explicit/relative-transformation.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/util/serialization.hh>
#include <pinocchio/serialization/se3.hpp>

#include "../src/explicit/input-configurations.hh"

namespace hpp {
namespace constraints {
namespace explicit_ {
LiegroupSpacePtr_t RelativePose::SE3(LiegroupSpace::SE3());
LiegroupSpacePtr_t RelativePose::R3xSO3(LiegroupSpace::R3xSO3());

RelativePosePtr_t RelativePose::create(
    const std::string& name, const DevicePtr_t& robot,
    const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
    const Transform3s& frame1, const Transform3s& frame2,
    ComparisonTypes_t comp, std::vector<bool> mask) {
  if (mask.empty()) {
    mask = std::vector<bool>(6, true);
  }
  RelativePose* ptr(new RelativePose(name, robot, joint1, joint2, frame1,
                                     frame2, comp, mask));
  RelativePosePtr_t shPtr(ptr);
  RelativePoseWkPtr_t wkPtr(shPtr);
  ptr->init(wkPtr);
  return shPtr;
}

RelativePosePtr_t RelativePose::createCopy(const RelativePosePtr_t& other) {
  RelativePose* ptr(new RelativePose(*other));
  RelativePosePtr_t shPtr(ptr);
  RelativePoseWkPtr_t wkPtr(shPtr);
  ptr->init(wkPtr);
  return shPtr;
}

ImplicitPtr_t RelativePose::copy() const { return createCopy(weak_.lock()); }

void RelativePose::outputValue(LiegroupElementRef result, vectorIn_t qin,
                               LiegroupElementConstRef rhs) const {
  vector_t rhsExplicit(explicitFunction()->outputSpace()->nv());
  implicitToExplicitRhs(rhs, rhsExplicit);
  explicitFunction()->value(result, qin);
  result += rhsExplicit;
}

void RelativePose::jacobianOutputValue(vectorIn_t qin,
                                       LiegroupElementConstRef f_value,
                                       LiegroupElementConstRef rhs,
                                       matrixOut_t jacobian) const {
  // Todo: suppress dynamic allocation.
  vector_t rhsExplicit(explicitFunction()->outputSpace()->nv());
  explicitFunction()->jacobian(jacobian, qin);
  if (rhs != rhs.space()->neutral()) {
    implicitToExplicitRhs(rhs, rhsExplicit);
    explicitFunction()
        ->outputSpace()
        ->dIntegrate_dq<pinocchio::DerivativeTimesInput>(f_value, rhsExplicit,
                                                         jacobian);
  }
}

void RelativePose::implicitToExplicitRhs(LiegroupElementConstRef implicitRhs,
                                         vectorOut_t explicitRhs) const {
  assert(*(implicitRhs.space()) == *R3xSO3 || *(implicitRhs.space()) == *SE3);
  assert(explicitRhs.size() == 6);
  // convert implicitRhs to Transform3s M1
  Transform3s M1(
      (Quaternion_t(implicitRhs.vector().tail<4>())).toRotationMatrix(),
      implicitRhs.vector().head<3>());
  // M2 = F_{2/J_2} M1 F_{2/J_2}^{-1}
  Transform3s M2(frame2_ * M1 * frame2_.inverse());
  // convert M2 to LiegroupElement p2
  vector7_t v2;
  v2.head<3>() = M2.translation();
  v2.tail<4>() = Quaternion_t(M2.rotation()).coeffs();
  LiegroupElement p2(v2, SE3);
  explicitRhs = p2 - SE3->neutral();
}

void RelativePose::explicitToImplicitRhs(vectorIn_t explicitRhs,
                                         LiegroupElementRef implicitRhs) const {
  assert(*(implicitRhs.space()) == *(R3xSO3) || *(implicitRhs.space()) == *SE3);
  assert(explicitRhs.size() == 6);
  // p1 = exp_{SE(3)} (explicitRhs)
  LiegroupElement p1(SE3->exp(explicitRhs));
  // convert p1 to Transform3s M1
  Transform3s M1((Quaternion_t(p1.vector().tail<4>())).toRotationMatrix(),
                 p1.vector().head<3>());
  // M2 = F_{2/J_2}^{-1} M1 F_{2/J_2}
  Transform3s M2(frame2_.inverse() * M1 * frame2_);
  // convert M2 to LiegroupElement p2
  vector7_t v2;
  implicitRhs.vector().head<3>() = M2.translation();
  implicitRhs.vector().tail<4>() = Quaternion_t(M2.rotation()).coeffs();
}

RelativePose::RelativePose(const std::string& name, const DevicePtr_t& robot,
                           const JointConstPtr_t& joint1,
                           const JointConstPtr_t& joint2,
                           const Transform3s& frame1, const Transform3s& frame2,
                           ComparisonTypes_t comp, std::vector<bool> mask)
    : Explicit(RelativeTransformationR3xSO3::create(name, robot, joint1, joint2,
                                                    frame1, frame2,
                                                    std::vector<bool>(6, true)),
               RelativeTransformation::create(name, robot, joint1, joint2,
                                              frame1, frame2),
               relativePose::inputConfVariables(robot, joint1, joint2),
               relativePose::jointConfInterval(joint2),
               relativePose::inputVelocityVariables(robot, joint1, joint2),
               relativePose::jointVelInterval(joint2), comp, mask),
      joint1_(joint1),
      joint2_(joint2),
      frame1_(frame1),
      frame2_(frame2) {
  assert(
      (*joint2->configurationSpace() == *pinocchio::LiegroupSpace::SE3()) ||
      (*joint2->configurationSpace() == *pinocchio::LiegroupSpace::R3xSO3()));
}

RelativePose::RelativePose(const RelativePose& other)
    : Explicit(other),
      joint1_(other.joint1_),
      joint2_(other.joint2_),
      frame1_(other.frame1_),
      frame2_(other.frame2_) {}

void RelativePose::init(RelativePoseWkPtr_t weak) {
  Explicit::init(weak);
  weak_ = weak;
}

template <class Archive>
void RelativePose::serialize(Archive& ar, const unsigned int version) {
  (void)version;
  ar& boost::serialization::make_nvp(
      "base", boost::serialization::base_object<Explicit>(*this));
  ar& BOOST_SERIALIZATION_NVP(frame1_);
  ar& BOOST_SERIALIZATION_NVP(frame2_);
  ar& BOOST_SERIALIZATION_NVP(weak_);
}

HPP_SERIALIZATION_IMPLEMENT(RelativePose);
}  // namespace explicit_
}  // namespace constraints
}  // namespace hpp

BOOST_CLASS_EXPORT_IMPLEMENT(hpp::constraints::explicit_::RelativePose)
