// Copyright (c) 2020, Airbus SAS and CNRS
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

#include <hpp/constraints/affine-function.hh>
#include <hpp/constraints/convex-shape-contact.hh>
#include <hpp/constraints/explicit/convex-shape-contact.hh>
#include <hpp/constraints/explicit/relative-pose.hh>
#include <hpp/pinocchio/liegroup-space.hh>

#include "../src/explicit/input-configurations.hh"

namespace hpp {
namespace constraints {
namespace explicit_ {
// Check that object contact surface
//   - belong to the same object (joint),
//   - are on a freeflyer joint.
// Throw in case these assertions are not true
static void checkContactSurfaces(const JointAndShapes_t& objectSurfaces) {
  JointPtr_t joint(0x0);
  for (auto js : objectSurfaces) {
    if (!joint) {
      joint = js.first;
      // Check that joint is a freeflyer
      if ((*joint->configurationSpace() != *pinocchio::LiegroupSpace::SE3()) &&
          (*joint->configurationSpace() !=
           *pinocchio::LiegroupSpace::R3xSO3())) {
        std::ostringstream os;
        os << "You are trying to build an explicit contact constraint but"
              " the joint that holds at least on object contact surface is "
              " not a freeflyer joint: "
           << joint->name();
        throw std::logic_error(os.str().c_str());
        if (joint->parentJoint()) {
          os << "You are trying to build an explicit contact constraint "
                "but the joint that holds at least on object contact "
                "surface: "
             << joint->name()
             << " is attached to another joint. This is not supported";
          throw std::logic_error(os.str().c_str());
        }
      }
    } else if (js.first != joint) {
      std::ostringstream os;
      os << "You are trying to build an explicit contact constraint "
            "but several joints hold object contact surface: "
         << joint->name() << " and " << js.first->name();
      throw std::logic_error(os.str().c_str());
    }
  }
}

ConvexShapeContactPtr_t ConvexShapeContact::create(
    const std::string& name, DevicePtr_t robot,
    const JointAndShapes_t& floorSurfaces,
    const JointAndShapes_t& objectSurfaces, const value_type& margin) {
  checkContactSurfaces(objectSurfaces);
  ConvexShapeContact* ptr(new ConvexShapeContact(name, robot, floorSurfaces,
                                                 objectSurfaces, margin));
  ConvexShapeContactPtr_t shPtr(ptr);
  ptr->init(shPtr);
  return shPtr;
}

ConvexShapeContact::Constraints_t
ConvexShapeContact::createConstraintAndComplement(
    const std::string& name, DevicePtr_t robot,
    const JointAndShapes_t& floorSurfaces,
    const JointAndShapes_t& objectSurfaces, const value_type& margin) {
  Constraints_t result;
  std::pair<hpp::constraints::ConvexShapeContactPtr_t,
            hpp::constraints::ConvexShapeContactComplementPtr_t>
      functions(ConvexShapeContactComplement::createPair(
          name, robot, floorSurfaces, objectSurfaces));
  functions.first->setNormalMargin(margin);
  // Contact constraint (= 0)
  std::get<0>(result) = Implicit::create(functions.first, 5 * EqualToZero);
  // Contact constraint complement (= rhs)
  std::get<1>(result) = Implicit::create(
      functions.second,
      ComparisonTypes_t(functions.second->outputSize(), constraints::Equality));
  std::get<2>(result) =
      create(name + "/hold", robot, floorSurfaces, objectSurfaces, margin);
  return result;
}

ConvexShapeContactPtr_t ConvexShapeContact::createCopy(
    const ConvexShapeContactPtr_t& other) {
  ConvexShapeContact* ptr(new ConvexShapeContact(*other));
  ConvexShapeContactPtr_t shPtr(ptr);
  ConvexShapeContactWkPtr_t wkPtr(shPtr);
  ptr->init(wkPtr);
  return shPtr;
}

ImplicitPtr_t ConvexShapeContact::copy() const {
  return createCopy(weak_.lock());
}

void ConvexShapeContact::outputValue(LiegroupElementRef result, vectorIn_t qin,
                                     LiegroupElementConstRef rhs) const {
  assert(HPP_DYNAMIC_PTR_CAST(ConvexShapeContactHold, functionPtr()));
  ConvexShapeContactHoldPtr_t f(
      HPP_STATIC_PTR_CAST(ConvexShapeContactHold, functionPtr()));
  std::size_t ifloor, iobject;
  LiegroupElement relativePoseRhs(LiegroupSpace::R3xSO3());
  f->complement()->computeRelativePoseRightHandSide(rhs, ifloor, iobject,
                                                    relativePoseRhs);
  // Extract input configuration of relative pose from qin
  Eigen::RowBlockIndices inputIndices(inputConf());
  vector_t q(f->inputSize());
  q.fill(sqrt(-1));
  inputIndices.lview(q) = qin;
  RelativePosePtr_t relativePose(pose_[ifloor * nFloor_ + iobject]);
  Eigen::RowBlockIndices relPosInputIndices(relativePose->inputConf());
  vector_t qinRelPose = relPosInputIndices.rview(q);
  assert(!qinRelPose.hasNaN());
  relativePose->outputValue(result, qinRelPose, relativePoseRhs);
}

void ConvexShapeContact::jacobianOutputValue(vectorIn_t qin,
                                             LiegroupElementConstRef,
                                             LiegroupElementConstRef rhs,
                                             matrixOut_t jacobian) const {
  assert(HPP_DYNAMIC_PTR_CAST(ConvexShapeContactHold, functionPtr()));
  ConvexShapeContactHoldPtr_t f(
      HPP_STATIC_PTR_CAST(ConvexShapeContactHold, functionPtr()));
  std::size_t ifloor, iobject;
  LiegroupElement relativePoseRhs(LiegroupSpace::R3xSO3());
  f->complement()->computeRelativePoseRightHandSide(rhs, ifloor, iobject,
                                                    relativePoseRhs);
  // Extract input configuration of relative pose from qin
  Eigen::RowBlockIndices inputIndices(inputConf());
  vector_t q(f->inputSize());
  q.fill(sqrt(-1));
  inputIndices.lview(q) = qin;
  RelativePosePtr_t relativePose(pose_[ifloor * nFloor_ + iobject]);
  Eigen::RowBlockIndices relPosInputIndices(relativePose->inputConf());
  vector_t qinRelPose = relPosInputIndices.rview(q);
  assert(!qinRelPose.hasNaN());
  LiegroupElement outputRelPose(LiegroupSpace::R3xSO3());
  relativePose->outputValue(outputRelPose, qinRelPose, relativePoseRhs);
  relativePose->jacobianOutputValue(qinRelPose, outputRelPose, relativePoseRhs,
                                    jacobian);
}

ConvexShapeContact::ConvexShapeContact(const std::string& name,
                                       DevicePtr_t robot,
                                       const JointAndShapes_t& floorSurfaces,
                                       const JointAndShapes_t& objectSurfaces,
                                       const value_type& margin)
    : Explicit(
          ConvexShapeContactHold::create(name, robot, floorSurfaces,
                                         objectSurfaces),
          ConstantFunction::create(
              pinocchio::LiegroupSpace::SE3()->neutral(),
              contact::inputSize(robot, floorSurfaces, objectSurfaces),
              contact::inputDerivSize(robot, floorSurfaces, objectSurfaces),
              name),
          contact::inputConfVariables(robot, floorSurfaces, objectSurfaces),
          relativePose::jointConfInterval(objectSurfaces.front().first),
          contact::inputVelocityVariables(robot, floorSurfaces, objectSurfaces),
          relativePose::jointVelInterval(objectSurfaces.front().first),
          (5 * EqualToZero << 3 * Equality), std::vector<bool>(8, true)) {
  assert(HPP_DYNAMIC_PTR_CAST(ConvexShapeContactHold, functionPtr()));
  ConvexShapeContactHoldPtr_t f(
      HPP_STATIC_PTR_CAST(ConvexShapeContactHold, functionPtr()));
  f->contactConstraint()->setNormalMargin(margin);
  const ConvexShapes_t& fs(f->contactConstraint()->floorContactSurfaces());
  const ConvexShapes_t& os(f->contactConstraint()->objectContactSurfaces());
  nFloor_ = fs.size();
  // Compute explicit relative poses
  for (std::size_t j = 0; j < fs.size(); ++j) {
    // move floor surface along x to take into account margin.
    Transform3f posInJoint(fs[j].positionInJoint());
    hppDout(info, "posInJoint" << posInJoint);
    posInJoint.translation() -= margin * posInJoint.rotation().col(0);
    hppDout(info, "posInJoint" << posInJoint);
    for (std::size_t i = 0; i < os.size(); ++i) {
      // Create explicit relative pose for each combination
      // (floor surface, object surface)
      pose_.push_back(
          RelativePose::create("", robot, fs[j].joint_, os[i].joint_,
                               posInJoint, os[i].positionInJoint(),
                               EqualToZero << 3 * Equality << 2 * EqualToZero,
                               std::vector<bool>(6, true)));
    }
  }
}
void ConvexShapeContact::init(ConvexShapeContactWkPtr_t weak) {
  Explicit::init(weak);
  weak_ = weak;
}
}  // namespace explicit_
}  // namespace constraints
}  // namespace hpp
