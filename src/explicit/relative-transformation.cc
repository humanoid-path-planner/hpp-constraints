// Copyright (c) 2017 - 2018, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
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

#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/weak_ptr.hpp>
#include <hpp/constraints/explicit/relative-transformation.hh>
#include <hpp/constraints/serialization.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/util/serialization.hh>
#include <pinocchio/serialization/se3.hpp>
#include <pinocchio/spatial/explog.hpp>
#include <pinocchio/spatial/skew.hpp>

#include "../serialization.hh"

namespace hpp {
namespace constraints {
namespace explicit_ {
namespace {
typedef JointJacobian_t::ConstNRowsBlockXpr<3>::Type ConstHalfJacobian_t;
inline ConstHalfJacobian_t omega(const JointJacobian_t& j) {
  return j.bottomRows<3>();
}
inline ConstHalfJacobian_t trans(const JointJacobian_t& j) {
  return j.topRows<3>();
}

void inputVariable(JointConstPtr_t joint, std::vector<bool>& conf,
                   std::vector<bool>& vel) {
  while (joint && joint->index() != 0) {
    for (size_type i = 0; i < joint->configSize(); ++i)
      conf[joint->rankInConfiguration() + i] =
          !conf[joint->rankInConfiguration() + i];
    for (size_type i = 0; i < joint->numberDof(); ++i)
      vel[joint->rankInVelocity() + i] = !vel[joint->rankInVelocity() + i];
    hppDout(info, "Adding joint " << joint->name() << " as input variable.");
    joint = joint->parentJoint();
  }
}

BlockIndex::segments_t vectorOfBoolToIntervals(std::vector<bool>& v) {
  BlockIndex::segments_t ret;
  for (std::size_t i = 0; i < v.size(); ++i)
    if (v[i]) ret.push_back(BlockIndex::segment_t(i, 1));
  BlockIndex::shrink(ret);
  return ret;
}

BlockIndex::segments_t jointConfInterval(JointConstPtr_t j) {
  return BlockIndex::segments_t(
      1, BlockIndex::segment_t(j->rankInConfiguration(), j->configSize()));
}
BlockIndex::segments_t jointVelInterval(JointConstPtr_t j) {
  return BlockIndex::segments_t(
      1, BlockIndex::segment_t(j->rankInVelocity(), j->numberDof()));
}
}  // namespace

RelativeTransformationPtr_t RelativeTransformation::create(
    const std::string& name, const DevicePtr_t& robot,
    const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
    const Transform3f& frame1, const Transform3f& frame2) {
  std::vector<bool> conf(robot->configSize(), false);
  std::vector<bool> vel(robot->numberDof(), false);
  inputVariable(joint1, conf, vel);
  inputVariable(joint2->parentJoint(), conf, vel);

  RelativeTransformation* ptr = new RelativeTransformation(
      name, robot, joint1, joint2, frame1, frame2,
      vectorOfBoolToIntervals(conf), jointConfInterval(joint2),
      vectorOfBoolToIntervals(vel), jointVelInterval(joint2));
  RelativeTransformationPtr_t shPtr(ptr);
  ptr->init(shPtr);
  return shPtr;
}

RelativeTransformation::RelativeTransformation(
    const std::string& name, const DevicePtr_t& robot,
    const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
    const Transform3f& frame1, const Transform3f& frame2,
    const segments_t inConf, const segments_t outConf, const segments_t inVel,
    const segments_t outVel, std::vector<bool> /*mask*/)
    : DifferentiableFunction(BlockIndex::cardinal(inConf),
                             BlockIndex::cardinal(inVel),
                             pinocchio::LiegroupSpace::SE3(), name),
      robot_(robot),
      parentJoint_(joint2->parentJoint()),
      joint1_(joint1),
      joint2_(joint2),
      frame1_(frame1),
      frame2_(frame2),
      inConf_(inConf),
      inVel_(inVel),
      outConf_(outConf),
      outVel_(outVel),
      F1inJ1_invF2inJ2_(frame1 * frame2.inverse()) {}

void RelativeTransformation::forwardKinematics(vectorIn_t arg) const {
  qsmall_ = inConf_.rview(robot_->currentConfiguration());
  if (qsmall_ != arg) {
    q_ = robot_->currentConfiguration();
    inConf_.lview(q_) = arg;
    robot_->currentConfiguration(q_);
  }
  robot_->computeForwardKinematics(pinocchio::JOINT_POSITION);
}

void RelativeTransformation::impl_compute(LiegroupElementRef result,
                                          vectorIn_t argument) const {
  forwardKinematics(argument);

  bool hasParent = (parentJoint_ && parentJoint_->index() > 0);

  // J1 * M1/J1 = J2 * M2/J2
  // J2 = J1 * M1/J1 * M2/J2^{-1}
  // J2 = J2_{parent} * T
  // T = J2_{parent}^{-1} * J2
  // T = J2_{parent}^{-1} * J1 * F1/J1 * F2/J2^{-1}
  Transform3f freeflyerPose;
  if (!joint1_)
    freeflyerPose = F1inJ1_invF2inJ2_;
  else
    freeflyerPose = joint1_->currentTransformation() * F1inJ1_invF2inJ2_;

  if (hasParent)
    freeflyerPose = parentJoint_->currentTransformation().actInv(freeflyerPose);

  freeflyerPose = joint2_->positionInParentFrame().actInv(freeflyerPose);

  typedef Transform3f::Quaternion Q_t;
  result.vector().head<3>() = freeflyerPose.translation();
  result.vector().tail<4>() = Q_t(freeflyerPose.rotation()).coeffs();
}

void RelativeTransformation::impl_jacobian(matrixOut_t jacobian,
                                           vectorIn_t arg) const {
  LiegroupElement result(outputSpace());
  impl_compute(result, arg);
  Configuration_t q(robot_->currentConfiguration());
  outConf_.lview(q) = result.vector();
  robot_->currentConfiguration(q);
  robot_->computeForwardKinematics(pinocchio::JOINT_POSITION | pinocchio::JACOBIAN);

  bool absolute = !joint1_;
  bool hasParent = (parentJoint_ && parentJoint_->index() > 0);

  static const JointJacobian_t Jabs;
  const JointJacobian_t& J1(absolute ? Jabs : joint1_->jacobian());
  // const JointJacobian_t& J2_parent (parentJoint_->jacobian());

  const matrix3_t& R1(absolute ? matrix3_t::Identity().eval()
                               : joint1_->currentTransformation().rotation());
  const matrix3_t& R2(joint2_->currentTransformation().rotation());
  const matrix3_t& R2_inParentFrame(
      joint2_->positionInParentFrame().rotation());

  const vector3_t& t1(absolute
                          ? vector3_t::Zero().eval()
                          : joint1_->currentTransformation().translation());

  matrix3_t cross1 = ::pinocchio::skew(
                (R1 * F1inJ1_invF2inJ2_.translation()).eval()),
            cross2;
  if (hasParent) {
    const vector3_t& t2_parent(
        parentJoint_->currentTransformation().translation());
    cross2 = ::pinocchio::skew((t2_parent - t1).eval());

    if (absolute)
      J2_parent_minus_J1_.noalias() = parentJoint_->jacobian();
    else
      J2_parent_minus_J1_.noalias() = parentJoint_->jacobian() - J1;
  } else {
    cross2 = -::pinocchio::skew(t1);
    // J2_parent_minus_J1_ = - J1;
  }

  // Express velocity of J1 * M1/J1 * M2/J2^{-1} in J2_{parent}.
  if (hasParent) {
    const matrix3_t& R2_parent(
        parentJoint_->currentTransformation().rotation());
    const JointJacobian_t& J2_parent(parentJoint_->jacobian());

    tmpJac_.noalias() =
        (R2_inParentFrame.transpose() * R2_parent.transpose()) *
        (cross1 * (omega(J2_parent_minus_J1_)) - cross2 * omega(J2_parent) -
         trans(J2_parent_minus_J1_));
    jacobian.topRows<3>() = inVel_.rview(tmpJac_);
  } else {
    if (absolute)
      jacobian.topRows<3>().setZero();
    else {
      tmpJac_.noalias() =
          R2.transpose() * ((-cross1 * R1) * omega(J1) + R1 * trans(J1));
      jacobian.topRows<3>() = inVel_.rview(tmpJac_);
    }
  }

  if (hasParent) {
    const matrix3_t& R2_parent(
        parentJoint_->currentTransformation().rotation());
    const JointJacobian_t& J2_parent(parentJoint_->jacobian());

    // J = p2RT2 * 0RTp2 * [ p2
    tmpJac_.noalias() = (R2.transpose() * R2_parent) * omega(J2_parent);
    if (!absolute) tmpJac_.noalias() -= (R2.transpose() * R1) * omega(J1);
    jacobian.bottomRows<3>() = inVel_.rview(tmpJac_);
  } else {
    if (absolute)
      jacobian.bottomRows<3>().setZero();
    else {
      tmpJac_.noalias() = (R2.transpose() * R1) * omega(J1);
      jacobian.bottomRows<3>() = inVel_.rview(tmpJac_);
    }
  }
}

template <class Archive>
void RelativeTransformation::serialize(Archive& ar,
                                       const unsigned int version) {
  (void)version;
  ar& boost::serialization::make_nvp(
      "base", boost::serialization::base_object<DifferentiableFunction>(*this));
  ar& BOOST_SERIALIZATION_NVP(robot_);
  internal::serialize_joint(ar, "parentJoint_", parentJoint_);
  internal::serialize_joint(ar, "joint1_", joint1_);
  internal::serialize_joint(ar, "joint2_", joint2_);
  ar& BOOST_SERIALIZATION_NVP(frame1_);
  ar& BOOST_SERIALIZATION_NVP(frame2_);

  ar& BOOST_SERIALIZATION_NVP(inConf_);
  ar& BOOST_SERIALIZATION_NVP(inVel_);
  ar& BOOST_SERIALIZATION_NVP(outConf_);
  ar& BOOST_SERIALIZATION_NVP(outVel_);
  ar& BOOST_SERIALIZATION_NVP(F1inJ1_invF2inJ2_);
  ar& BOOST_SERIALIZATION_NVP(weak_);
}

HPP_SERIALIZATION_IMPLEMENT(RelativeTransformation);
}  // namespace explicit_
}  // namespace constraints
}  // namespace hpp

BOOST_CLASS_EXPORT_IMPLEMENT(
    hpp::constraints::explicit_::RelativeTransformation)
