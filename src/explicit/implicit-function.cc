// Copyright (c) 2015 - 2018 LAAS-CNRS
// Authors: Florent Lamiraux, Joseph Mirabel
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

#include <Eigen/Geometry>
#include <hpp/constraints/differentiable-function.hh>
#include <hpp/constraints/explicit/implicit-function.hh>
#include <hpp/constraints/matrix-view.hh>
#include <hpp/constraints/serialization.hh>
#include <hpp/constraints/tools.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/liegroup-space.hh>
#include <hpp/pinocchio/serialization.hh>
#include <hpp/util/serialization.hh>

namespace hpp {
namespace constraints {
namespace explicit_ {
typedef Eigen::Map<Quaternion_t> QuatMap_t;
typedef Eigen::Map<const Quaternion_t> QuatConstMap_t;
typedef hpp::pinocchio::liegroup::VectorSpaceOperation<3, false> R3;
typedef hpp::pinocchio::liegroup::SpecialOrthogonalOperation<3> SO3;
typedef hpp::pinocchio::liegroup::CartesianProductOperation<R3, SO3> R3xSO3;
typedef pinocchio::liegroup::SpecialEuclideanOperation<3> SE3;
typedef hpp::pinocchio::LiegroupType LiegroupType;

struct JacobianVisitor : public boost::static_visitor<> {
  JacobianVisitor(vectorIn_t qOut, vectorIn_t f_qIn, matrixIn_t Jf,
                  const Eigen::MatrixBlocks<false, false>& outJacobian,
                  const Eigen::MatrixBlocks<false, false>& inJacobian,
                  matrixOut_t result)
      : qOut_(qOut),
        f_qIn_(f_qIn),
        Jf_(Jf),
        outJacobian_(outJacobian),
        inJacobian_(inJacobian),
        result_(result) {}

  template <typename LgT>
  void operator()(const LgT&);

  vectorIn_t qOut_, f_qIn_;
  matrixIn_t Jf_;
  const Eigen::MatrixBlocks<false, false>& outJacobian_;
  const Eigen::MatrixBlocks<false, false>& inJacobian_;
  matrixOut_t result_;
};  // struct JacobianVisitor

template <>
inline void JacobianVisitor::operator()<R3xSO3>(const R3xSO3&) {
  using Eigen::BlockIndex;
  using Eigen::MatrixBlocks;
  typedef hpp::constraints::BlockIndex BlockIndex;
  hppDout(info, "result_ = " << std::endl << result_);
  // Fill R^3 part
  assert(outJacobian_.nbRows() == 6);
  Eigen::MatrixBlocks<false, false> tmp(outJacobian_.block(0, 0, 3, 3));
  tmp.lview(result_) = matrix3_t::Identity();
  hppDout(info, "result_ = " << std::endl << result_);
  // extract 3 top rows of inJacobian_
  segments_t cols(inJacobian_.cols());
  MatrixBlocks<false, false> inJacobian(
      inJacobian_.block(0, 0, 3, BlockIndex::cardinal(cols)));
  inJacobian.lview(result_) = -Jf_.topRows<3>();
  hppDout(info, "result_ = " << std::endl << result_);
  // Fill SO(3) part
  // extract 3 bottom rows of inJacobian_
  inJacobian = inJacobian_.block(3, 0, 3, BlockIndex::cardinal(cols));
  // extract 3x3 bottom left part of outJacobian_
  MatrixBlocks<false, false> outJacobian(outJacobian_.block(3, 3, 3, 3));
  assert(qOut_.size() == 7);
  assert(f_qIn_.size() == 7);
  matrix3_t R_out(Eigen::Quaterniond(qOut_.tail<4>()).toRotationMatrix());
  matrix3_t R_f(Eigen::Quaterniond(f_qIn_.tail<4>()).toRotationMatrix());
  // \f$R_f^T R_{out}\f$
  matrix3_t R_f_T_R_out(R_f.transpose() * R_out);
  matrix3_t Jlog_R_f_T_R_out;
  vector3_t r;
  value_type theta;
  constraints::logSO3(R_f_T_R_out, theta, r);
  constraints::JlogSO3(theta, r, Jlog_R_f_T_R_out);
  outJacobian.lview(result_) = Jlog_R_f_T_R_out;
  hppDout(info, "result_ = " << std::endl << result_);
  inJacobian.lview(result_) =
      -Jlog_R_f_T_R_out * R_out.transpose() * R_f * Jf_.bottomRows<3>();
  hppDout(info, "result_ = " << std::endl << result_);
}

template <>
inline void JacobianVisitor::operator()<SE3>(const SE3&) {
  assert(outJacobian_.nbRows() == 6);
  // extract 3 top rows of inJacobian_
  assert(qOut_.size() == 7);
  assert(f_qIn_.size() == 7);
  matrix3_t Rout(Eigen::Quaterniond(qOut_.tail<4>()).toRotationMatrix());
  vector3_t pOut(qOut_.head<3>());
  Transform3s Mout(Rout, pOut);
  matrix3_t Rf(Eigen::Quaterniond(f_qIn_.tail<4>()).toRotationMatrix());
  vector3_t pf(f_qIn_.head<3>());
  Transform3s Mf(Rf, pf);
  // \f$Mf^{-1} M_{out}\f$
  Transform3s Mf_inverse_Mout(Mf.inverse() * Mout);
  matrix6_t Jlog_Mf_inverse_Mout;
  constraints::JlogSE3(Mf_inverse_Mout, Jlog_Mf_inverse_Mout);
  outJacobian_.lview(result_) = Jlog_Mf_inverse_Mout;
  JointJacobian_t inJ(6, Jf_.cols());
  inJ.topRows<3>().noalias() =
      (Rout.transpose() * ::pinocchio::skew(pOut - pf) * Rf) *
      Jf_.bottomRows<3>();
  inJ.topRows<3>().noalias() -= (Rout.transpose() * Rf) * Jf_.topRows<3>();
  inJ.bottomRows<3>().noalias() =
      (-Rout.transpose() * Rf) * Jf_.bottomRows<3>();
  inJacobian_.lview(result_) = Jlog_Mf_inverse_Mout * inJ;
}

template <typename LgT>
void JacobianVisitor::operator()(const LgT&) {
  outJacobian_.lview(result_).setIdentity();
  inJacobian_.lview(result_) = -Jf_;
}

ImplicitFunctionPtr_t ImplicitFunction::create(
    const LiegroupSpacePtr_t& configSpace,
    const DifferentiableFunctionPtr_t& function, const segments_t& inputConf,
    const segments_t& outputConf, const segments_t& inputVelocity,
    const segments_t& outputVelocity) {
  ImplicitFunction* ptr =
      new ImplicitFunction(configSpace, function, inputConf, outputConf,
                           inputVelocity, outputVelocity);
  return ImplicitFunctionPtr_t(ptr);
}

const DifferentiableFunctionPtr_t& ImplicitFunction::inputToOutput() const {
  return inputToOutput_;
}

ImplicitFunction::ImplicitFunction(const LiegroupSpacePtr_t& configSpace,
                                   const DifferentiableFunctionPtr_t& function,
                                   const segments_t& inputConf,
                                   const segments_t& outputConf,
                                   const segments_t& inputVelocity,
                                   const segments_t& outputVelocity)
    : DifferentiableFunction(configSpace->nq(), configSpace->nv(),
                             LiegroupSpace::Rn(function->outputSpace()->nv()),
                             "implicit " + function->name()),
      robot_(),
      inputToOutput_(function),
      inputConfIntervals_(inputConf),
      outputConfIntervals_(outputConf),
      inputDerivIntervals_(inputVelocity),
      outputDerivIntervals_(outputVelocity),
      outJacobian_(),
      inJacobian_(),
      f_qIn_(function->outputSpace()),
      qOut_(function->outputSpace()),
      result_(outputSpace()) {
  // Check input consistency
  // Each configuration variable is either input or output
  assert(function->inputSize() + function->outputSize() <= configSpace->nq());
  // Each velocity variable is either input or output
  assert(function->inputDerivativeSize() + function->outputDerivativeSize() <=
         configSpace->nv());
  qIn_.resize(function->inputSize());
  Jf_.resize(function->outputDerivativeSize(), function->inputDerivativeSize());
  assert(BlockIndex::cardinal(outputConf) == function->outputSize());
  // Sum of velocity output interval sizes equal function output
  // derivative size
  assert(BlockIndex::cardinal(outputVelocity) ==
         function->outputDerivativeSize());
  computeJacobianBlocks();

  activeParameters_.setConstant(false);
  activeDerivativeParameters_.setConstant(false);
  inputConfIntervals_.lview(activeParameters_.matrix()).setConstant(true);
  inputDerivIntervals_.lview(activeDerivativeParameters_.matrix())
      .setConstant(true);
  outputConfIntervals_.lview(activeParameters_.matrix()).setConstant(true);
  outputDerivIntervals_.lview(activeDerivativeParameters_.matrix())
      .setConstant(true);
}

void ImplicitFunction::impl_compute(LiegroupElementRef result,
                                    vectorIn_t argument) const {
  using Eigen::MatrixBlocks;
  hppDout(info, "argument=" << argument.transpose());
  // Store q_{output} in result
  qOut_.vector() = outputConfIntervals_.rview(argument);
  hppDout(info, "qOut_=" << qOut_);
  // fill in q_{input}
  qIn_ = inputConfIntervals_.rview(argument);
  hppDout(info, "qIn_=" << qIn_);
  // compute  f (q_{input}) -> output_
  inputToOutput_->value(f_qIn_, qIn_);
  hppDout(info, "f_qIn_=" << f_qIn_);
  result.vector() = qOut_ - f_qIn_;
  hppDout(info, "result=" << result);
}

void ImplicitFunction::impl_jacobian(matrixOut_t jacobian,
                                     vectorIn_t arg) const {
  jacobian.setZero();
  size_type iq = 0, iv = 0, nq, nv;
  std::size_t rank = 0;
  impl_compute(result_, arg);
  inputToOutput_->jacobian(Jf_, qIn_);
  hppDout(info, "Jf_=" << std::endl << Jf_);
  // Fill Jacobian by set of lines corresponding to the types of Lie group
  // that compose the outputspace of input to output function.
  for (std::vector<LiegroupType>::const_iterator it =
           inputToOutput_->outputSpace()->liegroupTypes().begin();
       it != inputToOutput_->outputSpace()->liegroupTypes().end(); ++it) {
    nq = inputToOutput_->outputSpace()->nq(rank);
    nv = inputToOutput_->outputSpace()->nv(rank);
    JacobianVisitor v(qOut_.vector().segment(iq, nq),
                      f_qIn_.vector().segment(iq, nq), Jf_.middleRows(iv, nv),
                      outJacobian_[rank], inJacobian_[rank],
                      jacobian.middleRows(iv, nv));
    boost::apply_visitor(v, *it);
    iq += nq;
    iv += nv;
    ++rank;
  }
}

void ImplicitFunction::computeJacobianBlocks() {
  segments_t remainingCols(outputDerivIntervals_.indices());
  segments_t cols;
  size_type iv = 0, nv;
  std::size_t rank = 0;
  for (std::vector<LiegroupType>::const_iterator it =
           inputToOutput_->outputSpace()->liegroupTypes().begin();
       it != inputToOutput_->outputSpace()->liegroupTypes().end(); ++it) {
    nv = inputToOutput_->outputSpace()->nv(rank);
    cols = BlockIndex::split(remainingCols, nv);
    segments_t rows(1, std::make_pair(iv, nv));
    outJacobian_.push_back(Eigen::MatrixBlocks<false, false>(rows, cols));
    inJacobian_.push_back(Eigen::MatrixBlocks<false, false>(
        rows, inputDerivIntervals_.indices()));
    //
    iv += nv;
    ++rank;
  }
}

template <class Archive>
void ImplicitFunction::serialize(Archive& ar, const unsigned int version) {
  (void)version;
  ar& boost::serialization::make_nvp(
      "base", boost::serialization::base_object<DifferentiableFunction>(*this));
  // ar & BOOST_SERIALIZATION_NVP(robot_);
  ar& BOOST_SERIALIZATION_NVP(inputToOutput_);
  ar& BOOST_SERIALIZATION_NVP(inputConfIntervals_);
  ar& BOOST_SERIALIZATION_NVP(outputConfIntervals_);
  ar& BOOST_SERIALIZATION_NVP(inputDerivIntervals_);
  ar& BOOST_SERIALIZATION_NVP(outputDerivIntervals_);
  ar& BOOST_SERIALIZATION_NVP(outJacobian_);
  ar& BOOST_SERIALIZATION_NVP(inJacobian_);

  if (!Archive::is_saving::value) {  // loading
    f_qIn_ = LiegroupElement(inputToOutput_->outputSpace());
    qOut_ = LiegroupElement(inputToOutput_->outputSpace());
    result_ = LiegroupElement(outputSpace());
    qIn_.resize(inputToOutput_->inputSize());
    Jf_.resize(inputToOutput_->outputDerivativeSize(),
               inputToOutput_->inputDerivativeSize());
    // should we recompute them instead of storing them ?
    // computeJacobianBlocks ();
  }
}

std::pair<JointConstPtr_t, JointConstPtr_t>
ImplicitFunction::dependsOnRelPoseBetween(DeviceConstPtr_t robot) const {
  assert(robot);
  // check that this is a constant function
  if (inputToOutput_->inputSize() != 0) {
    return std::pair<JointConstPtr_t, JointConstPtr_t>(nullptr, nullptr);
  }

  // get the joints involved in the output config
  JointConstPtr_t j1;
  // check that output interval matches with the config range of one joint
  if (outputConfIntervals_.rows().size() != 1) {
    return std::pair<JointConstPtr_t, JointConstPtr_t>(nullptr, nullptr);
  }
  segment_t row = outputConfIntervals_.rows()[0];
  j1 = robot->getJointAtConfigRank(row.first);
  if (!j1 || row.second != j1->configSize()) {
    return std::pair<JointConstPtr_t, JointConstPtr_t>(nullptr, nullptr);
  }

  JointConstPtr_t j2 = j1->parentJoint();
  size_type index1 = j1->index();
  size_type index2 = (j2 ? j2->index() : 0);
  if (index1 <= index2) {
    return std::pair<JointConstPtr_t, JointConstPtr_t>(j1, j2);
  } else {
    return std::pair<JointConstPtr_t, JointConstPtr_t>(j2, j1);
  }
}

HPP_SERIALIZATION_IMPLEMENT(ImplicitFunction);
}  // namespace explicit_
}  // namespace constraints
}  // namespace hpp

BOOST_CLASS_EXPORT_IMPLEMENT(hpp::constraints::explicit_::ImplicitFunction)
