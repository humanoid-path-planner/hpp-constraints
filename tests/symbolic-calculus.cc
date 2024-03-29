// Copyright (c) 2015, LAAS-CNRS
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

#define BOOST_TEST_MODULE SymbolicCalculus

#include <math.h>
#include <stdlib.h>

#include <Eigen/Geometry>
#include <boost/test/included/unit_test.hpp>
#include <hpp/constraints/affine-function.hh>
#include <hpp/constraints/symbolic-calculus.hh>
#include <hpp/constraints/symbolic-function.hh>
#include <limits>
#include <pinocchio/fwd.hpp>

using namespace hpp::constraints;

typedef SymbolicFunction<JointFrame> JointFrameFunction;

template <class ValueType = eigen::vector3_t,
          class JacobianType = JacobianMatrix>
class PointTesterT : public CalculusBase<PointTesterT<ValueType, JacobianType>,
                                         ValueType, JacobianType> {
 public:
  typedef CalculusBase<PointTesterT<ValueType, JacobianType>, ValueType,
                       JacobianType>
      Parent_t;
  struct DataWrapper {
    ValueType value;
    JacobianType jacobian;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };

  HPP_CONSTRAINTS_CB_CREATE1(PointTesterT, DataWrapper*)

  PointTesterT(DataWrapper* d = NULL) : datas(d) {}

  PointTesterT(const Parent_t& other)
      : Parent_t(other), datas(static_cast<const PointTesterT&>(other).datas) {}

  void impl_value(const ConfigurationIn_t) { this->value_ = datas->value; }
  void impl_jacobian(const ConfigurationIn_t) {
    this->jacobian_ = datas->jacobian;
  }

  DataWrapper* datas;
};

namespace crossProduct {
typedef Eigen::Matrix<value_type, 1, 6> Config;
typedef Eigen::Matrix<value_type, 3, 1> Value;
typedef Eigen::Matrix<value_type, 3, 6, Eigen::RowMajor> Jacobian;
typedef PointTesterT<Value, Jacobian> PointTester;
typedef PointTester::DataWrapper DataWrapper;

void setWrappers(Config cfg, DataWrapper* d1, DataWrapper* d2) {
  d1->value = cfg.leftCols(3);
  d2->value = cfg.rightCols(3);
  d1->jacobian = Jacobian::Zero();
  d2->jacobian = Jacobian::Zero();
  d1->jacobian.leftCols(3).setIdentity();
  d2->jacobian.rightCols(3).setIdentity();
}

void value(Config c, Value& v) {
  v[0] = c[1] * c[5] - c[2] * c[4];
  v[1] = -c[0] * c[5] + c[2] * c[3];
  v[2] = c[0] * c[4] - c[1] * c[3];
}

void jacobian(Config cfg, Eigen::Ref<Jacobian> j) {
  j(0, 0) = 0;
  j(0, 1) = cfg[5];
  j(0, 2) = -cfg[4];
  j(0, 3) = 0;
  j(0, 4) = -cfg[2];
  j(0, 5) = cfg[1];
  j(1, 0) = -cfg[5];
  j(1, 1) = 0;
  j(1, 2) = cfg[3];
  j(1, 3) = cfg[2];
  j(1, 4) = 0;
  j(1, 5) = -cfg[0];
  j(2, 0) = cfg[4];
  j(2, 1) = -cfg[3];
  j(2, 2) = 0;
  j(2, 3) = -cfg[1];
  j(2, 4) = cfg[0];
  j(2, 5) = 0;
}
}  // namespace crossProduct

BOOST_AUTO_TEST_CASE(CrossProductTest) {
  using namespace crossProduct;
  DataWrapper* d1 = new DataWrapper();
  DataWrapper* d2 = new DataWrapper();
  Traits<PointTester>::Ptr_t p1 = PointTester::create(d1),
                             p2 = PointTester::create(d2);
  typedef CrossProduct<PointTester, PointTester> CP_t;
  Traits<CP_t>::Ptr_t cp = p1 ^ p2;
  Value v;
  Jacobian j;
  vector_t unused;
  for (size_t i = 0; i < 100; i++) {
    Config cfg = Config::Random();
    cp->invalidate();
    setWrappers(cfg, d1, d2);
    cp->computeValue(unused);
    value(cfg, v);
    jacobian(cfg, j);
    BOOST_CHECK(cp->value().isApprox(v));
    cp->computeJacobian(unused);
    BOOST_CHECK(cp->jacobian().isApprox(j));
  }
  delete d1;
  delete d2;
}

namespace difference {
typedef Eigen::Matrix<value_type, 1, 6> Config;
typedef Eigen::Matrix<value_type, 3, 1> Value;
typedef Eigen::Matrix<value_type, 3, 6, Eigen::RowMajor> Jacobian;
typedef PointTesterT<Value, Jacobian> PointTester;
typedef PointTester::DataWrapper DataWrapper;

void setWrappers(Config cfg, DataWrapper* d1, DataWrapper* d2) {
  d1->value = cfg.leftCols(3);
  d2->value = cfg.rightCols(3);
  d1->jacobian = Jacobian::Identity();
  d2->jacobian = Jacobian::Zero();
  d2->jacobian.rightCols(3).setIdentity();
}

void value(Config c, Value& v) { v = c.leftCols(3) - c.rightCols(3); }

void jacobian(Config, Jacobian& j) {
  j.rightCols(3).setIdentity();
  j = -j;
  j.leftCols(3).setIdentity();
}
}  // namespace difference

BOOST_AUTO_TEST_CASE(DifferenceTest) {
  using namespace difference;
  DataWrapper* d1 = new DataWrapper();
  DataWrapper* d2 = new DataWrapper();
  Traits<PointTester>::Ptr_t p1 = PointTester::create(d1),
                             p2 = PointTester::create(d2);
  typedef Difference<PointTester, PointTester> D_t;
  Traits<D_t>::Ptr_t cp = p1 - p2;
  Value v;
  Jacobian j;
  vector_t unused;
  for (size_t i = 0; i < 100; i++) {
    Config cfg = Config::Random();
    cp->invalidate();
    setWrappers(cfg, d1, d2);
    cp->computeValue(unused);
    value(cfg, v);
    jacobian(cfg, j);
    BOOST_CHECK(cp->value().isApprox(v));
    cp->computeJacobian(unused);
    BOOST_CHECK(cp->jacobian().isApprox(j));
  }
  delete d1;
  delete d2;
}

namespace sum {
typedef Eigen::Matrix<value_type, 1, 6> Config;
typedef Eigen::Matrix<value_type, 3, 1> Value;
typedef Eigen::Matrix<value_type, 3, 6, Eigen::RowMajor> Jacobian;
typedef PointTesterT<Value, Jacobian> PointTester;
typedef PointTester::DataWrapper DataWrapper;

void setWrappers(Config cfg, DataWrapper* d1, DataWrapper* d2) {
  d1->value = cfg.leftCols(3);
  d2->value = cfg.rightCols(3);
  d1->jacobian = Jacobian::Identity();
  d2->jacobian = Jacobian::Zero();
  d2->jacobian.rightCols(3).setIdentity();
}

void value(Config c, Value& v) { v = c.leftCols(3) + c.rightCols(3); }

void jacobian(Config, Jacobian& j) {
  j.rightCols(3).setIdentity();
  j.leftCols(3).setIdentity();
}
}  // namespace sum

BOOST_AUTO_TEST_CASE(SumTest) {
  using namespace sum;
  DataWrapper* d1 = new DataWrapper();
  DataWrapper* d2 = new DataWrapper();
  Traits<PointTester>::Ptr_t p1 = PointTester::create(d1),
                             p2 = PointTester::create(d2);
  typedef Sum<PointTester, PointTester> D_t;
  Traits<D_t>::Ptr_t cp = p1 + p2;
  Value v;
  Jacobian j;
  vector_t unused;
  for (size_t i = 0; i < 100; i++) {
    Config cfg = Config::Random();
    cp->invalidate();
    setWrappers(cfg, d1, d2);
    cp->computeValue(unused);
    value(cfg, v);
    jacobian(cfg, j);
    BOOST_CHECK(cp->value().isApprox(v));
    cp->computeJacobian(unused);
    BOOST_CHECK(cp->jacobian().isApprox(j));
  }
  delete d1;
  delete d2;
}

namespace scalarMultiply {
typedef Eigen::Matrix<value_type, 1, 4> Config;
typedef Eigen::Matrix<value_type, 3, 1> Value;
typedef Eigen::Matrix<value_type, 3, 3, Eigen::RowMajor> Jacobian;
typedef PointTesterT<Value, Jacobian> PointTester;
typedef PointTester::DataWrapper DataWrapper;

void setWrappers(Config cfg, DataWrapper* d1, DataWrapper*) {
  d1->value = cfg.leftCols(3);
  d1->jacobian = Jacobian::Identity();
}

void value(Config c, Value& v) { v = c[3] * c.leftCols(3); }

void jacobian(Config c, Jacobian& j) {
  j.setIdentity();
  j = c[3] * j;
}
}  // namespace scalarMultiply

BOOST_AUTO_TEST_CASE(ScalarMultiplyTest) {
  using namespace scalarMultiply;
  DataWrapper* d1 = new DataWrapper();
  DataWrapper* d2 = new DataWrapper();
  Traits<PointTester>::Ptr_t p1 = PointTester::create(d1),
                             p2 = PointTester::create(d2);
  typedef ScalarMultiply<PointTester> D_t;
  value_type scalar = 3;
  Traits<D_t>::Ptr_t cp = scalar * p1;
  Value v;
  Jacobian j;
  vector_t unused;
  for (size_t i = 0; i < 100; i++) {
    Config cfg = Config::Random();
    cfg[3] = scalar;
    cp->invalidate();
    setWrappers(cfg, d1, d2);
    cp->computeValue(unused);
    value(cfg, v);
    jacobian(cfg, j);
    BOOST_CHECK(cp->value().isApprox(v));
    cp->computeJacobian(unused);
    BOOST_CHECK(cp->jacobian().isApprox(j));
  }
  delete d1;
  delete d2;
}

namespace scalarProduct {
typedef Eigen::Matrix<value_type, 1, 6> Config;
typedef Eigen::Matrix<value_type, 1, 1> OutputValue;
typedef Eigen::Matrix<value_type, 1, 6, Eigen::RowMajor> OutputJacobian;
typedef Eigen::Matrix<value_type, 3, 1> Value;
typedef Eigen::Matrix<value_type, 3, 6, Eigen::RowMajor> Jacobian;
typedef PointTesterT<Value, Jacobian> PointTester;
typedef PointTester::DataWrapper DataWrapper;

void setWrappers(Config cfg, DataWrapper* d1, DataWrapper* d2) {
  d1->value = cfg.leftCols(3);
  d2->value = cfg.rightCols(3);
  d1->jacobian = Jacobian::Zero();
  d2->jacobian = Jacobian::Zero();
  d1->jacobian.leftCols(3).setIdentity();
  d2->jacobian.rightCols(3).setIdentity();
}

void value(Config c, OutputValue& v) {
  v[0] = c.segment(0, 3).dot(c.segment(3, 3));
}

void jacobian(Config cfg, OutputJacobian& j) {
  j(0, 0) = cfg[3];
  j(0, 1) = cfg[4];
  j(0, 2) = cfg[5];
  j(0, 3) = cfg[0];
  j(0, 4) = cfg[1];
  j(0, 5) = cfg[2];
}
}  // namespace scalarProduct

BOOST_AUTO_TEST_CASE(ScalarProductTest) {
  using namespace scalarProduct;
  DataWrapper* d1 = new DataWrapper();
  DataWrapper* d2 = new DataWrapper();
  Traits<PointTester>::Ptr_t p1 = PointTester::create(d1),
                             p2 = PointTester::create(d2);
  typedef ScalarProduct<PointTester, PointTester> SP_t;
  Traits<SP_t>::Ptr_t sp = p1 * p2;
  OutputValue v;
  OutputJacobian j;
  vector_t unused;
  for (size_t i = 0; i < 100; i++) {
    Config cfg = Config::Random();
    sp->invalidate();
    setWrappers(cfg, d1, d2);
    sp->computeValue(unused);
    value(cfg, v);
    jacobian(cfg, j);
    BOOST_CHECK(v.isApprox(sp->value()));
    sp->computeJacobian(unused);
    BOOST_CHECK(sp->jacobian().isApprox(j));
  }
  delete d1;
  delete d2;
}

namespace matrixOfExp {
typedef Eigen::Matrix<value_type, 1, 6> Config;
typedef Eigen::Matrix<value_type, 3, 1> Value;
typedef Eigen::Matrix<value_type, 3, Eigen::Dynamic, Eigen::RowMajor> Jacobian;
typedef Eigen::Matrix<value_type, 6, 2> OutputValue;
typedef Eigen::Matrix<value_type, 6, 2 * 6, Eigen::RowMajor> OutputJacobian;
typedef Eigen::Matrix<value_type, 2, 1> RValue;
typedef Eigen::Matrix<value_type, 6, 6> OutputJacobianTimesRValue;
typedef PointTesterT<Value, Jacobian> PointTester;
typedef PointTester::DataWrapper DataWrapper;

void setWrappers(Config cfg, DataWrapper* d1, DataWrapper* d2) {
  d1->value = cfg.leftCols(3);
  d2->value = cfg.rightCols(3);
  d1->jacobian = Jacobian::Zero(3, 6);
  d2->jacobian = Jacobian::Zero(3, 6);
  d1->jacobian.leftCols(3).setIdentity();
  d2->jacobian.rightCols(3).setIdentity();
}

void value(Config c, OutputValue& v) {
  v.block(0, 0, 3, 1).transpose() = c.segment<3>(0);
  v.block(3, 0, 3, 1).transpose() = c.segment<3>(0).cross(c.segment<3>(3));
  v.block(0, 1, 3, 1).transpose() = c.segment(0, 3) - c.segment(3, 3);
  v.block(3, 1, 3, 1).transpose() = c.segment(0, 3) + c.segment(3, 3);
}

void jacobian(Config cfg, OutputJacobian& j) {
  j.block<3, 3>(0, 0).setIdentity();
  j.block<3, 3>(0, 3).setZero();

  crossProduct::jacobian(cfg, j.block<3, 6>(3, 0));
  // j.block (0,3,3,3).setZero();

  j.block<3, 3>(0, 6).setIdentity();
  j.block<3, 3>(0, 9).setIdentity();
  j.block<3, 3>(0, 9) *= -1;

  j.block<3, 3>(3, 6).setIdentity();
  j.block<3, 3>(3, 9).setIdentity();
}

void jacobianTimes(Eigen::Ref<const OutputJacobian> jin,
                   Eigen::Ref<const RValue> rvalue,
                   Eigen::Ref<OutputJacobianTimesRValue> jout) {
  jout.block(0, 0, 3, 6) =
      rvalue[0] * jin.block(0, 0, 3, 6) + rvalue[1] * jin.block(0, 6, 3, 6);
  jout.block(3, 0, 3, 6) =
      rvalue[0] * jin.block(3, 0, 3, 6) + rvalue[1] * jin.block(3, 6, 3, 6);
}
}  // namespace matrixOfExp

BOOST_AUTO_TEST_CASE(MatrixOfExpTest) {
  using namespace matrixOfExp;
  DataWrapper* d1 = new DataWrapper();
  DataWrapper* d2 = new DataWrapper();
  Traits<PointTester>::Ptr_t p1 = PointTester::create(d1),
                             p2 = PointTester::create(d2);
  typedef MatrixOfExpressions<Value, Jacobian> MoE_t;
  OutputValue v;
  OutputJacobian j;
  RValue rvalue;
  OutputJacobianTimesRValue jout, cache;
  MoE_t moe(v, j);
  moe.setSize(2, 2);
  moe.set(0, 0, p1);
  moe.set(1, 0, p1 ^ p2);
  moe.set(0, 1, p1 - p2);
  moe.set(1, 1, p1 + p2);
  vector_t unused;
  for (size_t i = 0; i < 100; i++) {
    Config cfg = Config::Random();
    moe.invalidate();
    setWrappers(cfg, d1, d2);
    moe.computeValue(unused);
    value(cfg, v);
    jacobian(cfg, j);
    BOOST_CHECK(moe.value().isApprox(v));
    moe.computeJacobian(unused);
    BOOST_CHECK(moe.jacobian().isApprox(j));
    for (size_t l = 0; l < 100; l++) {
      rvalue = RValue::Random();
      moe.jacobianTimes(unused, rvalue, cache);
      jacobianTimes(j, rvalue, jout);
      BOOST_CHECK(cache.isApprox(jout));
    }
  }
  delete d1;
  delete d2;
}

BOOST_AUTO_TEST_CASE(FunctionExpTest) {
  typedef Eigen::Matrix<value_type, 3, 1> Config;
  AffineFunctionPtr_t func(AffineFunction::create(matrix3_t::Identity()));
  hpp::shared_ptr<FunctionExp<AffineFunction> > funcExp =
      FunctionExp<AffineFunction>::create(func);
  for (size_t i = 0; i < 100; i++) {
    Config cfg = Config::Random();

    funcExp->invalidate();
    funcExp->computeValue(cfg);
    BOOST_CHECK(funcExp->value().isApprox(cfg));
    funcExp->computeJacobian(cfg);
    BOOST_CHECK(funcExp->jacobian().isIdentity());
  }
}
