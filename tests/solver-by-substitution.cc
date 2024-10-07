// Copyright (c) 2017, Joseph Mirabel
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

#define BOOST_TEST_MODULE SOLVER_BY_SUBSTITUTION
#include <../tests/util.hh>
#include <Eigen/Geometry>
#include <boost/test/unit_test.hpp>
#include <hpp/constraints/affine-function.hh>
#include <hpp/constraints/explicit/relative-pose.hh>
#include <hpp/constraints/generic-transformation.hh>
#include <hpp/constraints/solver/by-substitution.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/liegroup-element.hh>
#include <hpp/pinocchio/serialization.hh>
#include <hpp/pinocchio/simple-device.hh>
#include <pinocchio/algorithm/joint-configuration.hpp>
#include <sstream>

using hpp::constraints::AffineFunction;
using hpp::constraints::AffineFunctionPtr_t;
using hpp::constraints::ComparisonTypes_t;
using hpp::constraints::Configuration_t;
using hpp::constraints::ConstantFunction;
using hpp::constraints::ConstantFunctionPtr_t;
using hpp::constraints::DevicePtr_t;
using hpp::constraints::DifferentiableFunction;
using hpp::constraints::DifferentiableFunctionPtr_t;
using hpp::constraints::Equality;
using hpp::constraints::EqualToZero;
using hpp::constraints::Explicit;
using hpp::constraints::ExplicitConstraintSet;
using hpp::constraints::ExplicitPtr_t;
using hpp::constraints::Implicit;
using hpp::constraints::ImplicitPtr_t;
using hpp::constraints::JointPtr_t;
using hpp::constraints::LiegroupElement;
using hpp::constraints::LiegroupElementRef;
using hpp::constraints::LiegroupSpace;
using hpp::constraints::LockedJoint;
using hpp::constraints::matrix3_t;
using hpp::constraints::matrix_t;
using hpp::constraints::matrixOut_t;
using hpp::constraints::Orientation;
using hpp::constraints::RelativeTransformation;
using hpp::constraints::RelativeTransformationPtr_t;
using hpp::constraints::segment_t;
using hpp::constraints::segments_t;
using hpp::constraints::size_type;
using hpp::constraints::Transform3f;
using hpp::constraints::Transformation;
using hpp::constraints::value_type;
using hpp::constraints::vector3_t;
using hpp::constraints::vector_t;
using hpp::constraints::vectorIn_t;
using hpp::constraints::vectorOut_t;
using hpp::constraints::solver::BySubstitution;
using hpp::constraints::solver::lineSearch::Backtracking;
using hpp::constraints::solver::lineSearch::Constant;
using hpp::constraints::solver::lineSearch::ErrorNormBased;
using hpp::constraints::solver::lineSearch::FixedSequence;
using hpp::pinocchio::displayConfig;
using hpp::pinocchio::JACOBIAN;
using hpp::pinocchio::JOINT_POSITION;
using hpp::pinocchio::unittest::HumanoidRomeo;
using hpp::pinocchio::unittest::HumanoidSimple;
using hpp::pinocchio::unittest::makeDevice;
using hpp::pinocchio::unittest::ManipulatorArm2;

namespace saturation = hpp::constraints::solver::saturation;

matrix_t randomPositiveDefiniteMatrix(int N) {
  matrix_t A(matrix_t::Random(N, N));
  BOOST_REQUIRE((A.array() < 1).all());
  BOOST_REQUIRE((A.array() > -1).all());

  A = (A + A.transpose()) / 2;
  A += N * matrix_t::Identity(N, N);
  A /= N;
  return A;
}

const value_type test_precision = 1e-5;

//             x       y       z
template <int N1, int N2, int N3>
void test_quadratic() {
  const int N = N1 + N2 + N3;

  // (x y z) A (x y z)
  matrix_t A(randomPositiveDefiniteMatrix(N));
  Quadratic::Ptr_t quad(new Quadratic(A));

  // y = B * z
  matrix_t B(matrix_t::Random(N2, N3));
  const int Ninf = std::min(N2, N3);
  B.topLeftCorner(Ninf, Ninf) = randomPositiveDefiniteMatrix(Ninf);
  segments_t in;
  in.push_back(segment_t(N1 + N2, N3));
  segments_t out;
  out.push_back(segment_t(N1, N2));
  AffineFunctionPtr_t affine(AffineFunction::create(B));
  ExplicitPtr_t expl(
      Explicit::create(LiegroupSpace::Rn(N), affine, in, out, in, out));

  // Make solver
  BySubstitution solver(LiegroupSpace::Rn(N));
  solver.maxIterations(20);
  solver.errorThreshold(test_precision);
  solver.saturation(hpp::make_shared<saturation::Bounds>(-vector_t::Ones(N),
                                                         vector_t::Ones(N)));

  solver.add(Implicit::create(quad, ComparisonTypes_t(1, EqualToZero)));
  solver.add(expl);

  matrix_t M(N, N1 + N3);
  M << matrix_t::Identity(N1, N1), matrix_t::Zero(N1, N3),
      matrix_t::Zero(N2, N1), B, matrix_t::Zero(N3, N1),
      matrix_t::Identity(N3, N3);
  matrix_t Ar(M.transpose() * A * M);

  BOOST_CHECK_EQUAL(Ar.fullPivLu().rank(), N1 + N3);

  vector_t x(N);

  x.setZero();
  BOOST_CHECK(solver.isSatisfied(x));

  x.setRandom();
  SOLVER_CHECK_SOLVE(solver.solve<Backtracking>(x), SUCCESS);
  // EIGEN_VECTOR_IS_APPROX (x, vector_t::Zero(N));
  EIGEN_VECTOR_IS_APPROX(x.segment<N2>(N1), B * x.tail<N3>());
  BOOST_CHECK_SMALL(value_type(x.transpose() * A * x), test_precision);

  matrix_t expectedJ(1, N1 + N3), J(1, N1 + N3);

  x.setRandom();
  solver.explicitConstraintSet().solve(x);
  expectedJ =
      2 *
      solver.explicitConstraintSet().notOutArgs().rview(x).eval().transpose() *
      Ar;

  solver.computeValue<true>(x);
  solver.updateJacobian(x);
  solver.getReducedJacobian(J);

  EIGEN_IS_APPROX(expectedJ, J);
}

//             w       x       y       z
template <int N1, int N2, int N4, int N3>
void test_quadratic2() {
  const int N = N1 + N2 + N3 + N4;

  // (w x y z) A (w x y z)
  matrix_t A(randomPositiveDefiniteMatrix(N));
  Quadratic::Ptr_t quad(new Quadratic(A));

  // x = B * y
  matrix_t B(matrix_t::Random(N2, N3));
  const int Ninf = std::min(N2, N3);
  B.topLeftCorner(Ninf, Ninf) = randomPositiveDefiniteMatrix(Ninf);
  segments_t in1;
  in1.push_back(segment_t(N1 + N2, N3));
  segments_t out1;
  out1.push_back(segment_t(N1, N2));
  AffineFunctionPtr_t affine1(AffineFunction::create(B));
  ExplicitPtr_t expl1(
      Explicit::create(LiegroupSpace::Rn(N), affine1, in1, out1, in1, out1));

  // y = C * z
  matrix_t C(matrix_t::Random(N3, N4));
  const int Ninf2 = std::min(N3, N4);
  C.topLeftCorner(Ninf2, Ninf2) = randomPositiveDefiniteMatrix(Ninf2);
  segments_t in2;
  in2.push_back(segment_t(N1 + N2 + N3, N4));
  segments_t out2;
  out2.push_back(segment_t(N1 + N2, N3));
  AffineFunctionPtr_t affine2(AffineFunction::create(C));
  ExplicitPtr_t expl2(
      Explicit::create(LiegroupSpace::Rn(N), affine2, in2, out2, in2, out2));

  // Make solver
  BySubstitution solver(LiegroupSpace::Rn(N));
  solver.maxIterations(20);
  solver.errorThreshold(test_precision);
  solver.saturation(hpp::make_shared<saturation::Bounds>(-vector_t::Ones(N),
                                                         vector_t::Ones(N)));

  solver.add(Implicit::create(quad, ComparisonTypes_t(1, EqualToZero)));
  solver.add(expl1);
  solver.add(expl2);
  solver.explicitConstraintSetHasChanged();
  assert(solver.contains(expl1));
  assert(solver.contains(expl2));

  matrix_t M(N, N1 + N4);
  M << matrix_t::Identity(N1, N1), matrix_t::Zero(N1, N4),
      matrix_t::Zero(N2, N1), B * C, matrix_t::Zero(N3, N1), C,
      matrix_t::Zero(N4, N1), matrix_t::Identity(N4, N4);
  matrix_t Ar(M.transpose() * A * M);

  BOOST_CHECK_EQUAL(Ar.fullPivLu().rank(), N1 + N4);

  vector_t x(N);

  x.setZero();
  BOOST_CHECK(solver.isSatisfied(x));

  x.setRandom();
  SOLVER_CHECK_SOLVE(solver.solve<Backtracking>(x), SUCCESS);
  // SOLVER_CHECK_SOLVE (solver.solve<lineSearch::Constant>(x), SUCCESS);
  // EIGEN_VECTOR_IS_APPROX (x, vector_t::Zero(N));
  EIGEN_VECTOR_IS_APPROX(x.segment<N2>(N1), B * x.segment<N3>(N1 + N2));
  EIGEN_VECTOR_IS_APPROX(x.segment<N3>(N1 + N2),
                         C * x.segment<N4>(N1 + N2 + N3));
  BOOST_CHECK_SMALL(value_type(x.transpose() * A * x), test_precision);

  matrix_t expectedJ(1, N1 + N4), J(1, N1 + N4);

  x.setRandom();
  solver.explicitConstraintSet().solve(x);
  expectedJ =
      2 *
      solver.explicitConstraintSet().notOutArgs().rview(x).eval().transpose() *
      Ar;

  solver.computeValue<true>(x);
  solver.updateJacobian(x);
  solver.getReducedJacobian(J);

  EIGEN_IS_APPROX(expectedJ, J);
}

//             w       x       y       z
template <int N1, int N2, int N4, int N3>
void test_quadratic3() {
  const int N = N1 + N2 + N3 + N4;

  // x = B * (y z)
  matrix_t B(matrix_t::Random(N2, N3 + N4));
  const int Ninf = std::min(N2, N3 + N4);
  B.topLeftCorner(Ninf, Ninf) = randomPositiveDefiniteMatrix(Ninf);
  segments_t in1;
  in1.push_back(segment_t(N1 + N2, N3 + N4));
  segments_t out1;
  out1.push_back(segment_t(N1, N2));
  AffineFunctionPtr_t affine1(AffineFunction::create(B));
  ExplicitPtr_t expl1(
      Explicit::create(LiegroupSpace::Rn(N), affine1, in1, out1, in1, out1));

  // y = C * z
  matrix_t C(matrix_t::Random(N3, N4));
  const int Ninf2 = std::min(N3, N4);
  C.topLeftCorner(Ninf2, Ninf2) = randomPositiveDefiniteMatrix(Ninf2);
  segments_t in2;
  in2.push_back(segment_t(N1 + N2 + N3, N4));
  segments_t out2;
  out2.push_back(segment_t(N1 + N2, N3));
  AffineFunctionPtr_t affine2(AffineFunction::create(C));
  ExplicitPtr_t expl2(
      Explicit::create(LiegroupSpace::Rn(N), affine2, in2, out2, in2, out2));

  // z[0] = d
  vector_t d(vector_t::Random(1));
  segments_t in3;
  segments_t out3;
  out3.push_back(segment_t(N1 + N2 + N3, 1));
  ConstantFunctionPtr_t constant3(ConstantFunction::create(d, 0, 0));
  ExplicitPtr_t expl3(
      Explicit::create(LiegroupSpace::Rn(N), constant3, in3, out3, in3, out3));

  // (w x y z) A (w x y z)
  matrix_t A(randomPositiveDefiniteMatrix(N));
  Quadratic::Ptr_t quad(new Quadratic(A, -d[0]));

  // Make solver
  BySubstitution solver(LiegroupSpace::Rn(N));
  solver.maxIterations(20);
  solver.errorThreshold(test_precision);
  solver.saturation(hpp::make_shared<saturation::Bounds>(-vector_t::Ones(N),
                                                         vector_t::Ones(N)));

  solver.add(Implicit::create(quad, ComparisonTypes_t(1, Equality)));
  solver.add(expl1);
  solver.add(expl2);
  solver.add(expl3);
  BySubstitution copySolver(solver);

  matrix_t M(N, N1 + N4);
  M << matrix_t::Identity(N1, N1), matrix_t::Zero(N1, N4),
      matrix_t::Zero(N2, N1), B.leftCols(N3) * C + B.rightCols(N4),
      matrix_t::Zero(N3, N1), C, matrix_t::Zero(N4, N1),
      matrix_t::Identity(N4, N4);
  matrix_t P(N1 + N4, N1 + N4 - 1);
  P << matrix_t::Identity(N1, N1), matrix_t::Zero(N1, N4 - 1),
      matrix_t::Zero(1, N1 + N4 - 1), matrix_t::Zero(N4 - 1, N1),
      matrix_t::Identity(N4 - 1, N4 - 1);
  vector_t Xr_0(vector_t::Zero(N1 + N4));
  Xr_0[N1] = d[0];

  matrix_t Ar(M.transpose() * A * M);

  BOOST_CHECK_EQUAL(Ar.fullPivLu().rank(), N1 + N4);

  vector_t x(N);

  x.setRandom();
  SOLVER_CHECK_SOLVE(copySolver.solve<Backtracking>(x), SUCCESS);
  // SOLVER_CHECK_SOLVE (solver.solve<lineSearch::Constant>(x), SUCCESS);
  // EIGEN_VECTOR_IS_APPROX (x, vector_t::Zero(N));
  EIGEN_VECTOR_IS_APPROX(x.segment<N2>(N1), B * x.segment<N3 + N4>(N1 + N2));
  EIGEN_VECTOR_IS_APPROX(x.segment<N3>(N1 + N2),
                         C * x.segment<N4>(N1 + N2 + N3));
  BOOST_CHECK_SMALL(value_type(x.transpose() * A * x - d[0]), test_precision);

  matrix_t expectedJ(1, N1 + N4 - 1), J(1, N1 + N4 - 1);

  x.setRandom();
  copySolver.explicitConstraintSet().solve(x);
  expectedJ =
      2 *
      (P * copySolver.explicitConstraintSet().notOutArgs().rview(x).eval() +
       Xr_0)
          .transpose() *
      Ar * P;

  copySolver.computeValue<true>(x);
  copySolver.updateJacobian(x);
  copySolver.getReducedJacobian(J);

  EIGEN_IS_APPROX(expectedJ, J);
}

BOOST_AUTO_TEST_CASE(quadratic) {
  test_quadratic<3, 3, 3>();
  test_quadratic<5, 3, 4>();

  test_quadratic2<3, 3, 3, 3>();
  test_quadratic2<3, 4, 2, 6>();

  test_quadratic3<3, 3, 3, 3>();
  test_quadratic3<1, 4, 2, 6>();
}

void se3ToConfig(const Transform3f& oMi, vectorOut_t v) {
  assert(v.size() == 7);
  v.head<3>() = oMi.translation();
  Eigen::Map<Transform3f::Quaternion> q(v.tail<4>().data());
  q = oMi.rotation();
}

class Frame : public DifferentiableFunction {
 public:
  JointPtr_t joint_;

  Frame(JointPtr_t joint)
      : DifferentiableFunction(joint->robot()->configSize(),
                               joint->robot()->numberDof(),
                               LiegroupSpace::SE3(), "Frame"),
        joint_(joint) {}

  void impl_compute(LiegroupElementRef result, vectorIn_t arg) const {
    hpp::pinocchio::DeviceSync robot(joint_->robot());
    robot.currentConfiguration(arg);
    robot.computeForwardKinematics(JOINT_POSITION);

    const Transform3f& oMi = joint_->currentTransformation(robot.d());
    se3ToConfig(oMi, result.vector());
  }

  void impl_jacobian(matrixOut_t J, vectorIn_t arg) const {
    // finiteDifferenceCentral(J, arg, joint_->robot(), 1e-6);
    hpp::pinocchio::DeviceSync robot(joint_->robot());
    robot.currentConfiguration(arg);
    robot.computeForwardKinematics(JOINT_POSITION | JACOBIAN);

    J = joint_->jacobian(robot.d(), true);
  }
};

matrix3_t exponential(const vector3_t& aa) {
  matrix3_t R, xCross;
  xCross.setZero();
  xCross(1, 0) = +aa(2);
  xCross(0, 1) = -aa(2);
  xCross(2, 0) = -aa(1);
  xCross(0, 2) = +aa(1);
  xCross(2, 1) = +aa(0);
  xCross(1, 2) = -aa(0);
  R.setIdentity();
  value_type theta = aa.norm();
  if (theta < 1e-6) {
    R += xCross;
    R += 0.5 * xCross.transpose() * xCross;
  } else {
    R += sin(theta) / theta * xCross;
    R += 2 * std::pow(sin(theta / 2), 2) / std::pow(theta, 2) * xCross * xCross;
  }
  return R;
}

// Pose of a joint in root joint frame
//
// This differentiable function returns as ouput an element of R3xSO3
// corresponding to the pose of a joint J in the root joint frame of a robot.
// The input is the vector of configuration variables that move the joint with
// respect to the root joint.
class ExplicitTransformation : public DifferentiableFunction {
 public:
  JointPtr_t joint_;
  size_type in_, inDer_;
  // log_{SO(3)} (rootJoint^{-1}.joint)
  RelativeTransformationPtr_t rt_;

  // Constructor
  // joint: joint J the pose of which is computed with respect to root joint.
  // [in:in+l]         : interval of robot configuration variables that
  //                     modify the position of J with respect to the root
  //                     joint.
  // [inDer:inDer+lDer]: interval of robot velocity variables that modify the
  //                     position of J with respect to the root joint.
  ExplicitTransformation(JointPtr_t joint, size_type in, size_type l,
                         size_type inDer, size_type lDer)
      : DifferentiableFunction(l, lDer, LiegroupSpace::R3xSO3(),
                               "ExplicitTransformation"),
        joint_(joint),
        in_(in),
        inDer_(inDer) {
    rt_ = RelativeTransformation::create("RT", joint_->robot(),
                                         joint_->robot()->rootJoint(), joint_,
                                         Transform3f::Identity());
  }

  ExplicitConstraintSet::RowBlockIndices inArg() const {
    ExplicitConstraintSet::RowBlockIndices ret;
    ret.addRow(in_, inputSize());
    return ret;
  }

  ExplicitConstraintSet::RowBlockIndices outArg() const {
    ExplicitConstraintSet::RowBlockIndices ret;
    ret.addRow(0, 7);
    return ret;
  }

  ExplicitConstraintSet::ColBlockIndices inDer() const {
    ExplicitConstraintSet::ColBlockIndices ret;
    ret.addCol(inDer_, inputDerivativeSize());
    return ret;
  }

  ExplicitConstraintSet::RowBlockIndices outDer() const {
    ExplicitConstraintSet::RowBlockIndices ret;
    ret.addRow(0, 6);
    return ret;
  }
  // Fill input variables with arg. Other variables are set to neutral
  vector_t config(vectorIn_t arg) const {
    vector_t q = joint_->robot()->neutralConfiguration();
    q.segment(in_, inputSize()) = arg;
    return q;
    // joint_->robot()->currentConfiguration(q);
    // joint_->robot()->computeForwardKinematics();
  }

  // Compute relative position of joint_ wrt root joint
  // arg: vector of configuration variables that move joint_ in root joint.
  // result: R3xSO3 element containing the result
  void impl_compute(LiegroupElementRef result, vectorIn_t arg) const {
    // forwardKinematics(arg);
    LiegroupElement transform(LiegroupSpace::Rn(6));
    vector_t q = config(arg);
    rt_->value(transform, q);
    result.vector().head<3>() = transform.vector().head<3>();
    result.vector().tail<4>() =
        Eigen::Quaternion<value_type>(exponential(transform.vector().tail<3>()))
            .coeffs();

    // Transform3f tf1 = joint_->robot()->rootJoint()->currentTransformation();
    // Transform3f tf2 = joint_->currentTransformation();
    // Transform3f tf = tf2.inverse() * tf1;

    // result.head<3> = tf.translation();
    // result.tail<4> = Eigen::Quaternion<value_type>(tf.rotation());
  }

  void impl_jacobian(matrixOut_t jacobian, vectorIn_t arg) const {
    // forwardKinematics(arg);
    matrix_t J(6, rt_->inputDerivativeSize());
    vector_t q = config(arg);
    rt_->jacobian(J, q);

    inDer().rview(J).writeTo(jacobian);
  }
};

typedef hpp::shared_ptr<ExplicitTransformation> ExplicitTransformationPtr_t;

BOOST_AUTO_TEST_CASE(functions1) {
  BySubstitution solver(LiegroupSpace::R3());
  BySubstitution solver1(LiegroupSpace::R3());
  BySubstitution solver2(LiegroupSpace::R3());
  BySubstitution solver3(LiegroupSpace::R3());

  /// System:
  /// f (q1, q2) = 0
  /// h (    q2) = 0
  ///         q1 = g(q3)
  ///         q2 = C

  // f
  ImplicitPtr_t impl(Implicit::create(
      AffineFunction::create(matrix_t::Identity(2, 3)), 2 * EqualToZero));
  solver.add(impl);
  // Test inclusion of manifolds
  solver1.add(impl->copy());
  BOOST_CHECK(solver.definesSubmanifoldOf(solver));
  BOOST_CHECK(solver.definesSubmanifoldOf(solver1));
  BOOST_CHECK(solver1.definesSubmanifoldOf(solver));

  // q1 = g(q3)
  Eigen::Matrix<value_type, 1, 1> Jg;
  Jg(0, 0) = 1;
  Eigen::RowBlockIndices inArg;
  inArg.addRow(2, 1);
  Eigen::ColBlockIndices inDer;
  inDer.addCol(2, 1);
  Eigen::RowBlockIndices outArg;
  outArg.addRow(1, 1);
  segments_t in;
  in.push_back(segment_t(2, 1));
  segments_t out;
  out.push_back(segment_t(0, 1));
  AffineFunctionPtr_t affine(AffineFunction::create(matrix_t::Ones(1, 1)));
  ExplicitPtr_t expl(
      Explicit::create(LiegroupSpace::R3(), affine, in, out, in, out));
  solver.add(expl);
  // Test inclusion of manifolds
  BOOST_CHECK(solver.definesSubmanifoldOf(solver1));
  solver1.add(expl->copy());
  // q2 = C
  affine = AffineFunction::create(matrix_t(1, 0), vector_t::Zero(1));
  in.clear();
  out.clear();
  out.push_back(segment_t(1, 1));
  expl = Explicit::create(LiegroupSpace::R3(), affine, in, out, in, out);
  solver.add(expl);
  BOOST_CHECK_EQUAL(solver.reducedDimension(), 2);
  BOOST_CHECK(solver.definesSubmanifoldOf(solver1));

  // h
  matrix_t h(1, 3);
  h << 0, 1, 0;
  solver.add(Implicit::create(AffineFunction::create(h),
                              ComparisonTypes_t(1, EqualToZero)));
  BOOST_CHECK_EQUAL(solver.dimension(), 3);
  BOOST_CHECK_EQUAL(solver.reducedDimension(), 2);
  BOOST_CHECK(solver.definesSubmanifoldOf(solver1));

  segments_t impDof{segment_t(2, 1)};
  BOOST_CHECK_EQUAL(solver.implicitDof(), impDof);
}

BOOST_AUTO_TEST_CASE(functions2) {
  BySubstitution solver(LiegroupSpace::R3());

  /// System:
  /// f (q1, q3) = 0
  /// q2 = g(q3)
  Eigen::Matrix<value_type, 2, 3> Jf;
  Jf << 1, 0, 0, 0, 0, 1;
  solver.add(Implicit::create(AffineFunction::create(Jf), 2 * EqualToZero));

  Eigen::Matrix<value_type, 1, 1> Jg;
  Jg(0, 0) = 1;
  Eigen::RowBlockIndices inArg;
  inArg.addRow(2, 1);
  Eigen::ColBlockIndices inDer;
  inDer.addCol(2, 1);
  Eigen::RowBlockIndices outArg;
  outArg.addRow(1, 1);
  ExplicitPtr_t expl(Explicit::create(
      LiegroupSpace::R3(), AffineFunction::create(Jg), inArg.indices(),
      outArg.indices(), inDer.indices(), outArg.indices()));
  ImplicitPtr_t c1(expl);
  solver.add(expl);
  BOOST_CHECK_EQUAL(solver.dimension(), 2);

  // We add to the system h(q3) = 0
  /// f (q1, q3) = 0
  /// h (    q3) = 0
  /// q2 = g(q3)
  // This function should not be removed from the system.
  Eigen::Matrix<value_type, 1, 3> Jh;
  Jh << 0, 0, 1;
  ImplicitPtr_t impl(Implicit::create(AffineFunction::create(Jh),
                                      ComparisonTypes_t(1, EqualToZero)));
  ImplicitPtr_t c2(impl);
  c2->comparisonType(ComparisonTypes_t(1, Equality));
  solver.add(impl);
  BOOST_CHECK_EQUAL(solver.dimension(), 3);

  // We add to the system q3 = C
  // Function h should be removed, f should not.
  vector_t C(1);
  C(0) = 0;
  segments_t out;
  out.push_back(segment_t(2, 1));
  expl = Explicit::create(LiegroupSpace::R3(),
                          AffineFunction::create(matrix_t(1, 0), C),
                          segments_t(), out, segments_t(), out);
  ImplicitPtr_t c3(expl);
  c3->comparisonType(ComparisonTypes_t(1, Equality));
  solver.add(expl);

  BOOST_CHECK_EQUAL(solver.dimension(), 3);
  BOOST_CHECK_EQUAL(solver.reducedDimension(), 2);

  segments_t impDof = {segment_t(0, 1)};
  BOOST_CHECK_EQUAL(solver.implicitDof(), impDof);
  // test right hand side access by functions.
  vector_t rhs1(vector_t::Zero(c1->rightHandSideSize()));
  vector_t rhs2(vector_t::Random(c2->rightHandSideSize()));
  vector_t rhs3(vector_t::Random(c3->rightHandSideSize()));

  BOOST_CHECK(solver.rightHandSide(c1, rhs1));
  BOOST_CHECK(solver.rightHandSide(c2, rhs2));
  BOOST_CHECK(solver.rightHandSide(c3, rhs3));

  vector_t tmp1(c1->rightHandSideSize());
  tmp1.fill(sqrt(-1));
  vector_t tmp2(c2->rightHandSideSize());
  tmp2.fill(sqrt(-1));
  vector_t tmp3(c3->rightHandSideSize());
  tmp3.fill(sqrt(-1));

  BOOST_CHECK(solver.getRightHandSide(c1, tmp1));
  BOOST_CHECK(solver.getRightHandSide(c2, tmp2));
  BOOST_CHECK(solver.getRightHandSide(c3, tmp3));

  BOOST_CHECK(tmp1 == rhs1);
  BOOST_CHECK(tmp2 == rhs2);
  BOOST_CHECK(tmp3 == rhs3);
}

BOOST_AUTO_TEST_CASE(hybrid_solver) {
  DevicePtr_t device(makeDevice(HumanoidSimple));
  BOOST_REQUIRE(device);
  BOOST_CHECK_EQUAL(device->rootJoint()->positionInParentFrame(),
                    Transform3f::Identity());
  device->rootJoint()->lowerBound(0, -1);
  device->rootJoint()->lowerBound(1, -1);
  device->rootJoint()->lowerBound(2, -1);
  device->rootJoint()->upperBound(0, 1);
  device->rootJoint()->upperBound(1, 1);
  device->rootJoint()->upperBound(2, 1);
  JointPtr_t ee1 = device->getJointByName("rleg6_joint"),
             ee2 = device->getJointByName("lleg6_joint"),
             ee3 = device->getJointByName("larm6_joint");

  Configuration_t q0 = device->neutralConfiguration();
  device->currentConfiguration(q0);
  device->computeForwardKinematics(JOINT_POSITION);

  JointPtr_t lleg6Joint(device->getJointByName("lleg6_joint"));
  BOOST_CHECK_EQUAL(lleg6Joint->rankInConfiguration(), 12);
  BOOST_CHECK_EQUAL(lleg6Joint->rankInVelocity(), 11);
  BOOST_CHECK_EQUAL(lleg6Joint->configSize(), 1);
  BOOST_CHECK_EQUAL(lleg6Joint->numberDof(), 1);
  // Compute a configuration that satisfies the constaints.
  // Compute relative position of "lleg6_joint" wrt root
  Transform3f Mlleg6(lleg6Joint->currentTransformation());
  Transform3f Mroot(device->rootJoint()->currentTransformation());
  Transform3f M(Mroot.inverse() * Mlleg6);
  q0.segment<3>(0) = M.translation();
  q0.segment<4>(3) = Eigen::Quaternion<value_type>(M.rotation()).coeffs();

  BySubstitution solver(device->configSpace());
  solver.maxIterations(40);
  solver.errorThreshold(1e-4);
  solver.saturation(hpp::make_shared<saturation::Device>(device));

  device->currentConfiguration(q0);
  device->computeForwardKinematics(JOINT_POSITION);
  Transform3f tf1(ee1->currentTransformation());
  Transform3f tf2(ee2->currentTransformation());
  Transform3f tf3(ee3->currentTransformation());

  solver.add(Implicit::create(
      Orientation::create("Orientation lleg6_joint", device, ee2, tf2),
      3 * EqualToZero));
  solver.add(Implicit::create(
      Orientation::create("Orientation larm6_joint", device, ee3, tf3),
      3 * EqualToZero));

  BOOST_CHECK(solver.numberStacks() == 1);

  ExplicitTransformationPtr_t et;
  {
    et.reset(new ExplicitTransformation(lleg6Joint, 7, 6, 6, 6));
  }
  // Add an explicit constraint that computes the pose of the root joint (FF)
  // output variables are [0:7] for configurations and [0:6] for velocities
  // output value is equal to relative pose of "lleg6_joint" with respect to
  // the root joint.
  BOOST_CHECK(solver.explicitConstraintSet().add(Explicit::create(
                  device->configSpace(), et, et->inArg().indices(),
                  et->outArg().indices(), et->inDer().indices(),
                  et->outDer().indices())) >= 0);
  solver.explicitConstraintSetHasChanged();
  BOOST_TEST_MESSAGE(solver << '\n');

  BOOST_CHECK(solver.isSatisfied(q0));

  vector_t v(vector_t::Random(device->numberDof()));
  v *= .1;
  Configuration_t qrand(q0);
  LiegroupElement g(qrand, device->configSpace());
  g += v;
  qrand = g.vector();
  BOOST_CHECK_EQUAL(solver.solve<Backtracking>(qrand), BySubstitution::SUCCESS);
  qrand = g.vector();
  BOOST_CHECK_EQUAL(solver.solve<ErrorNormBased>(qrand),
                    BySubstitution::SUCCESS);
  qrand = g.vector();
  BOOST_CHECK_EQUAL(solver.solve<FixedSequence>(qrand),
                    BySubstitution::SUCCESS);
  BOOST_CHECK_EQUAL(solver.solve<Constant>(qrand), BySubstitution::SUCCESS);
}

BOOST_AUTO_TEST_CASE(by_substitution_serialization) {
  DevicePtr_t device(makeDevice(HumanoidSimple));
  BOOST_REQUIRE(device);
  device->rootJoint()->lowerBound(0, -1);
  device->rootJoint()->lowerBound(1, -1);
  device->rootJoint()->lowerBound(2, -1);
  device->rootJoint()->upperBound(0, 1);
  device->rootJoint()->upperBound(1, 1);
  device->rootJoint()->upperBound(2, 1);
  JointPtr_t ee1 = device->getJointByName("rleg5_joint"),
             ee2 = device->getJointByName("lleg5_joint"),
             ee3 = device->getJointByName("larm5_joint");

  Configuration_t q = device->currentConfiguration(),
                  qrand = ::pinocchio::randomConfiguration(device->model());

  BySubstitution solver(device->configSpace());
  solver.maxIterations(20);
  solver.errorThreshold(1e-3);
  solver.saturation(hpp::make_shared<saturation::Device>(device));

  device->currentConfiguration(q);
  device->computeForwardKinematics(JOINT_POSITION);
  Transform3f tf1(ee1->currentTransformation());
  Transform3f tf2(ee2->currentTransformation());
  Transform3f tf3(ee3->currentTransformation());

  solver.add(Implicit::create(
      Orientation::create("Orientation RAnkleRoll", device, ee2, tf2),
      3 * Equality));
  solver.add(Implicit::create(
      Orientation::create("Orientation LWristPitch", device, ee3, tf3),
      3 * Equality));
  solver.add(LockedJoint::create(ee1, ee1->configurationSpace()->neutral()));

  BOOST_CHECK(solver.numberStacks() == 1);

  std::stringstream ss;
  {
    hpp::serialization::xml_oarchive oa(ss);
    oa.insert(device->name(), device.get());
    oa << boost::serialization::make_nvp("solver", solver);
  }

  BOOST_TEST_MESSAGE(ss.str());

  BySubstitution r_solver(device->configSpace());
  {
    hpp::serialization::xml_iarchive ia(ss);
    ia.insert(device->name(), device.get());
    ia >> boost::serialization::make_nvp("solver", r_solver);
  }

  std::ostringstream ss_result, ss_expect;
  ss_expect << solver << '\n';
  ss_result << r_solver << '\n';
  BOOST_CHECK_EQUAL(ss_expect.str(), ss_result.str());
}

BOOST_AUTO_TEST_CASE(hybrid_solver_rhs) {
  using namespace hpp::constraints;

  DevicePtr_t device(makeDevice(HumanoidRomeo));
  BOOST_REQUIRE(device);

  Configuration_t q, qrand;

  JointPtr_t left = device->getJointByName("LWristPitch");
  TransformationR3xSO3::Ptr_t frame(TransformationR3xSO3::create(
      "LWristPitch", device, left, Transform3f::Identity()));
  Transformation::Ptr_t logFrame(Transformation::create(
      "LWristPitch", device, left, Transform3f::Identity()));

  // Check the logFrame if the log6 of frame.
  LiegroupElement valueFrame(frame->outputSpace()),
      logValFrame(logFrame->outputSpace());
  matrix_t Jframe(6, device->numberDof()), JlogFrame(6, device->numberDof()),
      expectedJlogFrame(6, device->numberDof());
  LiegroupElement neutral = frame->outputSpace()->neutral();
  for (int i = 0; i < 100; ++i) {
    q = ::pinocchio::randomConfiguration(device->model());

    frame->value(valueFrame, q);
    logFrame->value(logValFrame, q);

    vector_t expectedLog = hpp::pinocchio::log(valueFrame);

    EIGEN_VECTOR_IS_APPROX(expectedLog, logValFrame.vector());

    frame->jacobian(Jframe, q);
    logFrame->jacobian(expectedJlogFrame, q);

    JlogFrame = Jframe;
    frame->outputSpace()->dDifference_dq1<hpp::pinocchio::DerivativeTimesInput>(
        neutral.vector(), valueFrame.vector(), JlogFrame);

    EIGEN_IS_APPROX(expectedJlogFrame, JlogFrame);
  }

  // Check that the solver can handle constraints with R3xSO3 outputs.
  ImplicitPtr_t constraint(Implicit::create(frame, 6 * Equality));

  BySubstitution solver(device->configSpace());
  solver.maxIterations(20);
  solver.errorThreshold(1e-3);
  solver.saturation(hpp::make_shared<saturation::Device>(device));

  solver.add(constraint);

  BOOST_CHECK_EQUAL(solver.rightHandSideSize(), 7);

  for (int i = 0; i < 100; ++i) {
    q = ::pinocchio::randomConfiguration(device->model()),

    device->currentConfiguration(q);
    device->computeForwardKinematics(JOINT_POSITION);
    Transform3f tf_expected(left->currentTransformation());
    vector_t rhs_expected(7), rhs(7);
    se3ToConfig(tf_expected, rhs_expected);

    solver.rightHandSideFromConfig(q);
    rhs = solver.rightHandSide();
    SE3CONFIG_IS_APPROX(rhs_expected, rhs);
    solver.getRightHandSide(constraint, rhs);
    SE3CONFIG_IS_APPROX(rhs_expected, rhs);

    solver.rightHandSideFromConfig(constraint, q);
    rhs = solver.rightHandSide();
    SE3CONFIG_IS_APPROX(rhs_expected, rhs);

    solver.rightHandSide(rhs_expected);
    rhs = solver.rightHandSide();
    SE3CONFIG_IS_APPROX(rhs_expected, rhs);

    solver.rightHandSide(constraint, rhs_expected);
    rhs = solver.rightHandSide();
    SE3CONFIG_IS_APPROX(rhs_expected, rhs);

    BOOST_CHECK_EQUAL(solver.solve<FixedSequence>(q), BySubstitution::SUCCESS);

    BySubstitution::Status status;
    for (int j = 0; j < 100; ++j) {
      qrand = ::pinocchio::randomConfiguration(device->model());
      status = solver.solve<FixedSequence>(qrand);
      if (status == BySubstitution::SUCCESS) break;
    }
    BOOST_CHECK_EQUAL(status, BySubstitution::SUCCESS);

    if (status == BySubstitution::SUCCESS) {
      device->currentConfiguration(qrand);
      device->computeForwardKinematics(JOINT_POSITION);
      Transform3f tf_result(left->currentTransformation());

      Transform3f id = tf_expected.actInv(tf_result);
      BOOST_CHECK_MESSAGE(id.isIdentity(1e-3), "Right hand side is different:\n"
                                                   << tf_result << '\n'
                                                   << tf_expected << '\n'
                                                   << id);
    }
  }
}

BOOST_AUTO_TEST_CASE(rightHandSide) {
  for (size_type i = 0; i < 1000; ++i) {
    size_type N(10);
    matrix_t A(randomPositiveDefiniteMatrix((int)N));
    AffineFunctionPtr_t affine(AffineFunction::create(A));
    vector_t b(vector_t::Random(N));
    ComparisonTypes_t comp(N, Equality);
    comp[1] = comp[3] = comp[5] = comp[6] = EqualToZero;
    b[1] = b[3] = b[5] = b[6] = 0;

    BySubstitution solver(LiegroupSpace::Rn(N));
    solver.maxIterations(20);
    solver.errorThreshold(test_precision);
    // Create constraint with various comparison types
    ImplicitPtr_t constraint(Implicit::create(affine, comp));
    solver.add(constraint);
    solver.rightHandSide(constraint, b);
    vector_t b1(N);
    solver.getRightHandSide(constraint, b1);
    BOOST_CHECK(b == b1);
    // Check resolution
    vector_t x(vector_t::Random(N));
    solver.solve(x);
    vector_t error(A * x - b);
    BOOST_CHECK_MESSAGE(error.norm() < test_precision,
                        "Error threshold exceeded. Error is "
                            << error.transpose() << ", norm " << error.norm()
                            << ". Precision is " << test_precision);
  }
}

BOOST_AUTO_TEST_CASE(rightHandSideFromConfig) {
  // Create a kinematic chain
  DevicePtr_t device = hpp::pinocchio::unittest::makeDevice(HumanoidSimple);
  JointPtr_t root = device->rootJoint(),
             ee1 = device->getJointByName("lleg5_joint"),
             ee2 = device->getJointByName("rleg5_joint");
  BOOST_REQUIRE(device);

  ComparisonTypes_t comp1(EqualToZero << Equality << EqualToZero << Equality
                                      << EqualToZero << Equality);
  assert(comp1[0] == EqualToZero);
  assert(comp1[2] == EqualToZero);
  assert(comp1[4] == EqualToZero);
  ComparisonTypes_t comp2(2 * Equality << 2 * EqualToZero << 2 * Equality);
  assert(comp2[0] == Equality);
  assert(comp2[1] == Equality);
  assert(comp2[2] == EqualToZero);
  assert(comp2[3] == EqualToZero);
  assert(comp2[4] == Equality);
  assert(comp2[5] == Equality);
  // Create two relative transformation constraints
  Transform3f tf1(Transform3f::Identity());
  vector3_t u;
  u << 0, -.2, 0;
  Transform3f tf2(Transform3f::Identity());
  tf2.translation(u);

  DifferentiableFunctionPtr_t h(
      RelativeTransformation::create("RelativeTransformation", device, ee1, ee2,
                                     tf1, tf2, std::vector<bool>(6, true)));
  ImplicitPtr_t c1(Implicit::create(h, comp1));
  u << 1.2, 0, -1;
  tf2.translation(u);
  ImplicitPtr_t c2(hpp::constraints::explicit_::RelativePose::create(
      "Transformation", device, JointPtr_t(), root, tf2, tf1, comp2,
      std::vector<bool>(6, true)));

  BySubstitution solver(device->configSpace());
  solver.maxIterations(20);
  solver.errorThreshold(test_precision);
  solver.add(c1);
  solver.add(c2);
  //           0
  //           rhs [1]
  // f1 (q) =  0
  //           rhs [3]
  //           0
  //           rhs [5]
  //
  //           rhs [6]
  //           0
  // f2 (q) =  0
  //           rhs [9]
  //           rhs [10]
  //           rhs [11]
  //           rhs [12]
  for (size_type i = 0; i < 1000; ++i) {
    Configuration_t q = ::pinocchio::randomConfiguration(device->model());
    bool success;
    // Set right hand side for both constraints from random configuration
    success = solver.rightHandSideFromConfig(c1, q);
    BOOST_CHECK(success);
    success = solver.rightHandSideFromConfig(c2, q);
    BOOST_CHECK(success);
    // Store right hand side for each constraint
    vector_t rhs1(6);
    rhs1.setZero();
    vector_t rhs1_(6);
    rhs1_.setZero();
    vector_t rhs2(7);
    rhs2.setZero();
    vector_t rhs2_(7);
    rhs2_.setZero();
    success = solver.getRightHandSide(c1, rhs1);
    BOOST_CHECK(success);
    success = solver.getRightHandSide(c2, rhs2);
    BOOST_CHECK(success);
    // Set right hand side for both constraints from other random configuration
    q = ::pinocchio::randomConfiguration(device->model());
    success = solver.rightHandSideFromConfig(c1, q);
    BOOST_CHECK(success);
    success = solver.rightHandSideFromConfig(c2, q);
    BOOST_CHECK(success);
    // Set right hand side from stored values
    success = solver.rightHandSide(c1, rhs1);
    BOOST_CHECK(success);
    success = solver.rightHandSide(c2, rhs2);
    BOOST_CHECK(success);
    // Get right hand side for each constraint and compare to stored values
    success = solver.getRightHandSide(c1, rhs1_);
    BOOST_CHECK(success);
    success = solver.getRightHandSide(c2, rhs2_);
    BOOST_CHECK(success);
    BOOST_CHECK((rhs1 - rhs1_).norm() < 1e-10);
    BOOST_CHECK((rhs2 - rhs2_).norm() < 1e-10);
  }
}

BOOST_AUTO_TEST_CASE(merge) {
  // Create a kinematic chain
  DevicePtr_t device = hpp::pinocchio::unittest::makeDevice(HumanoidSimple);
  JointPtr_t root = device->rootJoint(),
             ee1 = device->getJointByName("lleg5_joint"),
             ee2 = device->getJointByName("rleg5_joint");
  BOOST_REQUIRE(device);
  Configuration_t q(device->configSpace()->neutral().vector());
  ComparisonTypes_t comp1(6 * Equality);
  comp1[0] = comp1[2] = comp1[4] = EqualToZero;
  ComparisonTypes_t comp2(6 * Equality);
  comp2[1] = comp2[2] = EqualToZero;
  // Create two relative transformation constraints
  Transform3f tf1(Transform3f::Identity());
  vector3_t u;
  u << 0, -.2, 0;
  Transform3f tf2(Transform3f::Identity());
  tf2.translation(u);
  DifferentiableFunctionPtr_t h(RelativeTransformation::create(
      "RelativeTransformation", device, ee1, ee2, tf1, tf2));
  ImplicitPtr_t c1(Implicit::create(h, comp1));
  u << 1.2, 0, -1;
  tf2.translation(u);
  ImplicitPtr_t c2(
      LockedJoint::create(ee1, ee1->configurationSpace()->neutral()));

  ImplicitPtr_t c3(hpp::constraints::explicit_::RelativePose::create(
      "Transformation root", device, JointPtr_t(), root, tf2, tf1, comp2));

  BySubstitution solver1(device->configSpace());
  BySubstitution solver2(device->configSpace());
  solver1.maxIterations(20);
  solver1.errorThreshold(test_precision);
  solver1.add(c1->copy());
  solver1.add(c2->copy());
  solver2.add(c1->copy());
  solver2.add(c3->copy());
  // copy and merge solvers
  BySubstitution solver3(solver1);
  BySubstitution solver4(solver2);

  solver3.merge(solver4);

  BOOST_CHECK(solver3.numericalConstraints().size() == 3);
  BOOST_CHECK(solver3.contains(c1));
  BOOST_CHECK(solver3.contains(c2));
  BOOST_CHECK(solver3.contains(c3));
  BOOST_CHECK(solver3.rightHandSideFromConfig(c1, q));
  BOOST_CHECK(solver3.rightHandSideFromConfig(c2, q));
  BOOST_CHECK(solver3.rightHandSideFromConfig(c3, q));

  // Check computation of errors by constraint
  // order is c1 c2 c3
  vector_t error(solver3.errorSize());
  vector_t error1(c1->function().outputSpace()->nv());
  vector_t error2(c2->function().outputSpace()->nv());
  vector_t error3(c3->function().outputSpace()->nv());
  for (int i = 0; i < 2; ++i) {
    bool satisfied(solver3.isSatisfied(q, error));
    bool found;
    bool satisfied1(solver3.isConstraintSatisfied(c1, q, error1, found));
    BOOST_CHECK(found);
    bool satisfied2(solver3.isConstraintSatisfied(c2, q, error2, found));
    BOOST_CHECK(found);
    bool satisfied3(solver3.isConstraintSatisfied(c3, q, error3, found));
    BOOST_CHECK(found);
    BOOST_CHECK_EQUAL(satisfied, (satisfied1 && satisfied2 && satisfied3));
    size_type row = 0, nRows = error1.size();
    BOOST_CHECK(error1 == error.segment(row, nRows));
    row += nRows;
    nRows = error2.size();
    BOOST_CHECK(error2 == error.segment(row, nRows));
    row += nRows;
    nRows = error3.size();
    BOOST_CHECK(error3 == error.segment(row, nRows));

    q = ::pinocchio::randomConfiguration(device->model());
  }

  BySubstitution solver5(device->configSpace());
  solver5.add(c3);
  BOOST_CHECK(solver5.contains(c3->copy()));
}
