// Copyright (c) 2015, Joseph Mirabel
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

#include "hpp/constraints/qp-static-stability.hh"

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <limits>

#include "hpp/constraints/tools.hh"

namespace hpp {
namespace constraints {
namespace {
std::size_t forceDatasToNbContacts(
    const std::vector<ConvexShapeContact::ForceData>& fds) {
  std::size_t nb = 0;
  for (std::vector<ConvexShapeContact::ForceData>::const_iterator it =
           fds.begin();
       it != fds.end(); ++it)
    nb += it->points.size();
  return nb;
}
}  // namespace

const Eigen::Matrix<value_type, 6, 1> QPStaticStability::Gravity =
    (Eigen::Matrix<value_type, 6, 1>() << 0, 0, -1, 0, 0, 0).finished();
const Eigen::Matrix<value_type, 6, 1> QPStaticStability::MinusGravity =
    (Eigen::Matrix<value_type, 6, 1>() << 0, 0, +1, 0, 0, 0).finished();

QPStaticStability::QPStaticStability(const std::string& name,
                                     const DevicePtr_t& robot,
                                     const Contacts_t& contacts,
                                     const CenterOfMassComputationPtr_t& com)
    : DifferentiableFunction(robot->configSize(), robot->numberDof(),
                             LiegroupSpace::R1(), name),
      Zeros(new qpOASES::real_t[contacts.size()]),
      nWSR(40),
      robot_(robot),
      nbContacts_(contacts.size()),
      com_(com),
      H_(nbContacts_, nbContacts_),
      G_(nbContacts_),
      qp_((qpOASES::int_t)nbContacts_, qpOASES::HST_SEMIDEF),
      phi_(Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero(6, nbContacts_),
           Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero(
               6, nbContacts_ * robot->numberDof())),
      primal_(vector_t::Zero(nbContacts_)),
      dual_(vector_t::Zero(nbContacts_)) {
  VectorMap_t zeros(Zeros, nbContacts_);
  zeros.setZero();

  qpOASES::Options options;
  qp_.setOptions(options);

  qp_.setPrintLevel(qpOASES::PL_NONE);
  phi_.setSize(2, nbContacts_);
  Traits<PointCom>::Ptr_t OG = PointCom::create(com);
  for (std::size_t i = 0; i < contacts.size(); ++i) {
    Traits<PointInJoint>::Ptr_t OP2 = PointInJoint::create(
        contacts[i].joint, contacts[i].point, robot->numberDof());
    Traits<VectorInJoint>::Ptr_t n2 = VectorInJoint::create(
        contacts[i].joint, contacts[i].normal, robot->numberDof());

    phi_(0, i) = n2;
    phi_(1, i) = (OG - OP2) ^ n2;
  }
}

QPStaticStability::QPStaticStability(const std::string& name,
                                     const DevicePtr_t& robot,
                                     const std::vector<ForceData>& contacts,
                                     const CenterOfMassComputationPtr_t& com)
    : DifferentiableFunction(robot->configSize(), robot->numberDof(), 1, name),
      Zeros(new qpOASES::real_t[forceDatasToNbContacts(contacts)]),
      nWSR(40),
      robot_(robot),
      nbContacts_(forceDatasToNbContacts(contacts)),
      com_(com),
      H_(nbContacts_, nbContacts_),
      G_(nbContacts_),
      qp_((qpOASES::int_t)nbContacts_, qpOASES::HST_SEMIDEF),
      phi_(Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero(6, nbContacts_),
           Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero(
               6, nbContacts_ * robot->numberDof())),
      primal_(vector_t::Zero(nbContacts_)),
      dual_(vector_t::Zero(nbContacts_)) {
  VectorMap_t zeros(Zeros, nbContacts_);
  zeros.setZero();

  qpOASES::Options options;
  qp_.setOptions(options);

  qp_.setPrintLevel(qpOASES::PL_NONE);
  phi_.setSize(2, nbContacts_);
  Traits<PointCom>::Ptr_t OG = PointCom::create(com);
  std::size_t col = 0;
  for (std::size_t i = 0; i < contacts.size(); ++i) {
    Traits<VectorInJoint>::Ptr_t n = VectorInJoint::create(
        contacts[i].joint, contacts[i].normal, robot->numberDof());
    for (std::size_t j = 0; j < contacts[i].points.size(); ++j) {
      Traits<PointInJoint>::Ptr_t OP = PointInJoint::create(
          contacts[i].joint, contacts[i].points[j], robot->numberDof());

      phi_(0, col) = n;
      phi_(1, col) = (OG - OP) ^ n;
      col++;
    }
  }
}

QPStaticStabilityPtr_t QPStaticStability::create(
    const std::string& name, const DevicePtr_t& robot,
    const Contacts_t& contacts, const CenterOfMassComputationPtr_t& com) {
  return QPStaticStabilityPtr_t(
      new QPStaticStability(name, robot, contacts, com));
}

QPStaticStabilityPtr_t QPStaticStability::create(
    const std::string& name, const DevicePtr_t& robot,
    const std::vector<ForceData>& contacts,
    const CenterOfMassComputationPtr_t& com) {
  return QPStaticStabilityPtr_t(
      new QPStaticStability(name, robot, contacts, com));
}

QPStaticStabilityPtr_t QPStaticStability::create(
    const DevicePtr_t& robot, const Contacts_t& contacts,
    const CenterOfMassComputationPtr_t& com) {
  return create("QPStaticStability", robot, contacts, com);
}

void QPStaticStability::impl_compute(LiegroupElementRef result,
                                     ConfigurationIn_t argument) const {
  robot_->currentConfiguration(argument);
  robot_->computeForwardKinematics(pinocchio::JOINT_POSITION);

  phi_.invalidate();
  phi_.computeValue(argument);
  // phi_.computeSVD (argument);

  qpOASES::returnValue ret = solveQP(result.vector());
  if (ret != qpOASES::SUCCESSFUL_RETURN) {
    hppDout(error, "QP could not be solved. Error is " << ret);
  }
  if (!checkQPSol()) {
    hppDout(error, "QP solution does not satisfies the constraints");
  }
}

void QPStaticStability::impl_jacobian(matrixOut_t jacobian,
                                      ConfigurationIn_t argument) const {
  robot_->currentConfiguration(argument);
  robot_->computeForwardKinematics(pinocchio::JOINT_POSITION |
                                   pinocchio::JACOBIAN);

  phi_.invalidate();
  // phi_.computeSVD (argument);
  phi_.computeJacobian(argument);

  vector_t res(1);
  qpOASES::returnValue ret = solveQP(res);
  if (ret != qpOASES::SUCCESSFUL_RETURN) {
    hppDout(error, "QP could not be solved. Error is " << ret);
  }
  if (!checkQPSol()) {
    hppDout(error, "QP solution does not satisfies the constraints");
  }
  if (!checkStrictComplementarity()) {
    hppDout(error,
            "Strict complementary slackness does not hold. "
            "Jacobian WILL be wrong.");
  }

  // TODO preallocate this.
  matrix_t JT_phi_F(nbContacts_, robot_->numberDof());
  matrix_t J_F(6, robot_->numberDof());
  phi_.jacobianTransposeTimes(argument, phi_.value() * primal_, JT_phi_F);
  phi_.jacobianTimes(argument, primal_, J_F);

  jacobian = 0.5 * primal_.transpose() * JT_phi_F +
             (0.5 * phi_.value() * primal_ + Gravity).transpose() * J_F;
}

inline qpOASES::returnValue QPStaticStability::solveQP(
    vectorOut_t result) const {
  // TODO: Use the SVD to solve a smaller quadratic problem
  // Try to find a positive solution
  using qpOASES::SUCCESSFUL_RETURN;

  H_ = phi_.value().transpose() * phi_.value();
  G_ = phi_.value().transpose() * Gravity;

  qpOASES::int_t nwsr = nWSR;
  qp_.reset();
  qp_.setHessianType(qpOASES::HST_SEMIDEF);
  qpOASES::returnValue ret;
  ret = qp_.init(H_.data(), G_.data(), Zeros, 0, nwsr, 0);
  qp_.getPrimalSolution(primal_.data());
  qp_.getDualSolution(dual_.data());
  result[0] = 2 * qp_.getObjVal() + MinusGravity.squaredNorm();
  return ret;
}

bool QPStaticStability::checkQPSol() const {
  return (primal_.array() >= -1e-8).all();
}

bool QPStaticStability::checkStrictComplementarity() const {
  qpOASES::real_t eps = qp_.getOptions().boundTolerance;
  return ((primal_.array() > eps && dual_.cwiseAbs().array() <= eps) ||
          (dual_.array() > eps && primal_.cwiseAbs().array() <= eps))
      .all();
}
}  // namespace constraints
}  // namespace hpp
