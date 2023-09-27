// Copyright (c) 2014, LAAS-CNRS
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

#include "hpp/constraints/static-stability.hh"

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/liegroup-element.hh>
#include <limits>

#include "hpp/constraints/tools.hh"

namespace hpp {
namespace constraints {

using hpp::pinocchio::LiegroupElement;

const value_type StaticStability::G = 9.81;
const Eigen::Matrix<value_type, 6, 1> StaticStability::Gravity =
    (Eigen::Matrix<value_type, 6, 1>() << 0, 0, -1, 0, 0, 0).finished();

StaticStability::StaticStability(const std::string& name,
                                 const DevicePtr_t& robot,
                                 const Contacts_t& contacts,
                                 const CenterOfMassComputationPtr_t& com)
    : DifferentiableFunction(robot->configSize(), robot->numberDof(),
                             contacts.size() + 6, name),
      robot_(robot),
      contacts_(contacts),
      com_(com),
      phi_(Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero(6,
                                                              contacts.size()),
           Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero(
               6, contacts.size() * robot->numberDof())),
      u_(contacts.size()),
      uMinus_(contacts.size()),
      v_(contacts.size()),
      uDot_(contacts.size(), robot->numberDof()),
      uMinusDot_(contacts.size(), robot->numberDof()),
      vDot_(contacts.size(), robot->numberDof()),
      lambdaDot_(robot->numberDof()) {
  phi_.setSize(2, contacts.size());
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

StaticStabilityPtr_t StaticStability::create(
    const std::string& name, const DevicePtr_t& robot,
    const Contacts_t& contacts, const CenterOfMassComputationPtr_t& com) {
  return StaticStabilityPtr_t(new StaticStability(name, robot, contacts, com));
}

StaticStabilityPtr_t StaticStability::create(
    const DevicePtr_t& robot, const Contacts_t& contacts,
    const CenterOfMassComputationPtr_t& com) {
  return create("StaticStability", robot, contacts, com);
}

void StaticStability::impl_compute(LiegroupElementRef result,
                                   ConfigurationIn_t argument) const {
  robot_->currentConfiguration(argument);
  robot_->computeForwardKinematics(pinocchio::JOINT_POSITION);

  phi_.invalidate();

  phi_.computeSVD(argument);

  const Eigen::Matrix<value_type, 6, 1> G = -1 * Gravity;
  u_.noalias() = phi_.svd().solve(G);

  if (computeUminusAndV(u_, uMinus_, v_)) {
    // value_type lambda, unused_lMax; size_type iMax, iMin;
    // findBoundIndex (u_, v_, lambda, &iMin, unused_lMax, &iMax);
    value_type lambda = 1;

    result.vector().segment(0, contacts_.size()) = u_ + lambda * v_;
  } else {
    result.vector().segment(0, contacts_.size()) = u_;
  }
  result.vector().segment<6>(contacts_.size()) = Gravity + phi_.value() * u_;
}

void StaticStability::impl_jacobian(matrixOut_t jacobian,
                                    ConfigurationIn_t argument) const {
  robot_->currentConfiguration(argument);
  robot_->computeForwardKinematics(pinocchio::JOINT_POSITION |
                                   pinocchio::JACOBIAN);

  phi_.invalidate();

  phi_.computeSVD(argument);
  phi_.computeJacobian(argument);
  phi_.computePseudoInverse(argument);

  const Eigen::Matrix<value_type, 6, 1> G = -1 * Gravity;
  u_.noalias() = phi_.svd().solve(G);
  phi_.computePseudoInverseJacobian(argument, G);
  uDot_.noalias() = phi_.pinvJacobian();

  jacobian.block(0, 0, contacts_.size(), robot_->numberDof()).noalias() = uDot_;

  if (computeUminusAndV(u_, uMinus_, v_)) {
    matrix_t S = -matrix_t::Identity(u_.size(), u_.size());
    S.diagonal() = 1 * (u_.array() >= 0).select(0, -vector_t::Ones(u_.size()));

    // value_type lambda, unused_lMax; size_type iMax, iMin;
    // findBoundIndex (u_, v_, lambda, &iMin, unused_lMax, &iMax);
    // if (lambda < Eigen::NumTraits <value_type>::dummy_precision())
    // return;
    value_type lambda = 1;

    computeVDot(argument, uMinus_, S.diagonal(), uDot_, uMinusDot_, vDot_);

    // computeLambdaDot (u_, v_, iMin, uDot_, vDot_, lambdaDot_);

    // jacobian.block (0, 0, contacts_.size(), robot_->numberDof()).noalias ()
    // += lambda * vDot_ + v_ * lambdaDot_.transpose ();

    jacobian.block(0, 0, contacts_.size(), robot_->numberDof()).noalias() +=
        lambda * vDot_;
  }

  phi_.jacobianTimes(
      argument, u_,
      jacobian.block(contacts_.size(), 0, 6, robot_->numberDof()));
  phi_.computePseudoInverseJacobian(argument, Gravity);
  jacobian.block(contacts_.size(), 0, 6, robot_->numberDof()) +=
      -phi_.value() * phi_.pinvJacobian();
}

void StaticStability::findBoundIndex(vectorIn_t u, vectorIn_t v,
                                     value_type& lambdaMin, size_type* iMin,
                                     value_type& lambdaMax, size_type* iMax) {
  // Be carefull when v has small values.
  // Consider them as 0
  const value_type eps = Eigen::NumTraits<value_type>::dummy_precision();
  vector_t lambdas =
      (v.array().cwiseAbs() < eps).select(0, -u.cwiseQuotient(v));
  lambdaMin = (v.array() > eps).select(lambdas, 0).maxCoeff(iMin);
  lambdaMax = (v.array() < eps).select(lambdas, 0).minCoeff(iMax);
}

bool StaticStability::computeUminusAndV(vectorIn_t u, vectorOut_t uMinus,
                                        vectorOut_t v) const {
  using namespace hpp::pinocchio;

  uMinus.noalias() = (u.array() >= 0).select(0, -u);

  if (uMinus.isZero()) return false;

  size_type rank = phi_.svd().rank();
  v.noalias() = getV2<MoE_t::SVD_t>(phi_.svd(), rank) *
                (getV2<MoE_t::SVD_t>(phi_.svd(), rank).adjoint() * uMinus);
  // v.noalias() = uMinus;
  // v.noalias() -= getV1 <MoE_t::SVD_t> (phi_.svd()) *
  // ( getV1 <MoE_t::SVD_t> (phi_.svd()).adjoint() * uMinus );
  return true;
}

void StaticStability::computeVDot(const ConfigurationIn_t arg,
                                  vectorIn_t uMinus, vectorIn_t S,
                                  matrixIn_t uDot, matrixOut_t uMinusDot,
                                  matrixOut_t vDot) const {
  using namespace hpp::pinocchio;

  size_type rank = phi_.svd().rank();
  uMinusDot.noalias() = S.asDiagonal() * uDot;
  vDot.noalias() = uMinusDot;
  vDot.noalias() -=
      getV1<MoE_t::SVD_t>(phi_.svd(), rank) *
      (getV1<MoE_t::SVD_t>(phi_.svd(), rank).adjoint() * uMinusDot);

  // TODO: preallocate this matrix
  Eigen::Matrix<value_type, 6, Eigen::Dynamic> JphiTimesUMinus(
      6, robot_->numberDof());
  phi_.jacobianTimes(arg, uMinus, JphiTimesUMinus);
  vDot.noalias() -= phi_.pinv() * JphiTimesUMinus;

  phi_.computePseudoInverseJacobian(arg, phi_.value() * uMinus);
  vDot.noalias() -= phi_.pinvJacobian();
}

void StaticStability::computeLambdaDot(vectorIn_t u, vectorIn_t v,
                                       const std::size_t i0, matrixIn_t uDot,
                                       matrixIn_t vDot,
                                       vectorOut_t lambdaDot) const {
  if (std::abs(v(i0)) < Eigen::NumTraits<value_type>::dummy_precision())
    lambdaDot.setZero();
  else {
    lambdaDot.noalias() = -uDot.row(i0) / v(i0);
    lambdaDot.noalias() += (u(i0) / (v(i0) * v(i0))) * vDot.row(i0);
  }
}
}  // namespace constraints
}  // namespace hpp
