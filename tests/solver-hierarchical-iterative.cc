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

#define BOOST_TEST_MODULE HIERARCHICAL_ITERATIVE_SOLVER
#include <../tests/util.hh>
#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>
#include <functional>
#include <hpp/constraints/affine-function.hh>
#include <hpp/constraints/generic-transformation.hh>
#include <hpp/constraints/implicit.hh>
#include <hpp/constraints/solver/hierarchical-iterative.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/liegroup-space.hh>
#include <hpp/pinocchio/simple-device.hh>
#include <pinocchio/algorithm/joint-configuration.hpp>

using namespace hpp::constraints;
namespace saturation = hpp::constraints::solver::saturation;

const value_type test_precision = 1e-5;

#define VECTOR2(x0, x1) ((hpp::constraints::vector_t(2) << x0, x1).finished())

using Eigen::VectorXi;
using hpp::pinocchio::LiegroupSpace;

template <typename LineSearch>
struct test_base {
  solver::HierarchicalIterative solver;
  LineSearch ls;

  test_base(const size_type& d) : solver(LiegroupSpace::Rn(d)) {
    solver.maxIterations(20);
    solver.errorThreshold(test_precision);
    solver.saturation(hpp::make_shared<saturation::Bounds>(vector_t::Zero(2),
                                                           vector_t::Ones(2)));
  }

  vector_t success(value_type x0, value_type x1) {
    vector_t x(VECTOR2(x0, x1));
    BOOST_CHECK_EQUAL(solver.solve(x, ls),
                      solver::HierarchicalIterative::SUCCESS);
    return x;
  }

  vector_t failure(value_type x0, value_type x1) {
    vector_t x(VECTOR2(x0, x1));
    BOOST_CHECK_PREDICATE(
        std::not_equal_to<solver::HierarchicalIterative::Status>(),
        (solver.solve(x, ls))(solver::HierarchicalIterative::SUCCESS));
    return x;
  }
};

template <typename LineSearch = solver::lineSearch::Constant>
struct test_quadratic : test_base<LineSearch> {
  test_quadratic(const matrix_t& A) : test_base<LineSearch>(A.cols()) {
    // Find (x, y)
    // s.t. a * x^2 + b * y^2 - 1 = 0
    //      0 <= x <= 1
    //      0 <= y <= 1
    BOOST_REQUIRE_EQUAL(A.rows(), A.cols());
    BOOST_TEST_MESSAGE(A);
    Quadratic::Ptr_t f(new Quadratic(A, -1));

    this->solver.add(
        Implicit::create(
            f, ComparisonTypes_t(f->outputDerivativeSize(), Equality)),
        0);
    BOOST_CHECK(this->solver.numberStacks() == 1);
  }
};

BOOST_AUTO_TEST_CASE(quadratic) {
  matrix_t A(2, 2);

  A << 1, 0, 0, 1;
  test_quadratic<> test(A);
  BOOST_CHECK_EQUAL(test.failure(0, 0), VECTOR2(0, 0));
  test.success(0.1, 0);
  test.success(0, 0.1);
  test.success(0.5, 0.5);

  A << 2, 0, 0, 2;
  test = test_quadratic<>(A);
  test.success(0.1, 0);
  test.success(0, 0.1);
  test.success(0.5, 0.5);

  A << 0.5, 0, 0, 0.5;
  test = test_quadratic<>(A);
  // This is exact because of the saturation
  BOOST_CHECK_EQUAL(test.success(1, 0.001),
                    VECTOR2(1, 1));  // Slide on the border x = 1
  BOOST_CHECK_EQUAL(test.success(0.001, 1),
                    VECTOR2(1, 1));  // Slide on the border y = 1

  A << 0.75, 0, 0, 0.75;
  test_quadratic<solver::lineSearch::FixedSequence> test4(A);
  // This is not exact because the solver does not saturate.
  EIGEN_VECTOR_IS_APPROX(
      test4.success(1, 0.1),
      VECTOR2(1., 1 / sqrt(3)));  // Slide on the border x = 1
  EIGEN_VECTOR_IS_APPROX(
      test4.success(0.1, 1),
      VECTOR2(1 / sqrt(3), 1.));  // Slide on the border y = 1
  // There is an overshoot. To overcome this, the Hessian of the function should
  // be obtained.
  EIGEN_VECTOR_IS_NOT_APPROX(
      test4.success(1, 0.001),
      VECTOR2(1., 1 / sqrt(3)));  // Slide on the border x = 1
  EIGEN_VECTOR_IS_NOT_APPROX(
      test4.success(0.001, 1),
      VECTOR2(1 / sqrt(3), 1.));  // Slide on the border y = 1

  // Ellipsoid: computations are approximative
  A << 0.5, 0, 0, 2;
  test_quadratic<solver::lineSearch::FixedSequence> test1(A);
  BOOST_CHECK_EQUAL(test1.success(1, 0.5),
                    VECTOR2(1., 0.5));  // Slide on the border x = 1
  EIGEN_VECTOR_IS_APPROX(test1.success(1, 0.1),
                         VECTOR2(1., 0.5));  // Slide on the border x = 1
  EIGEN_VECTOR_IS_APPROX(test1.success(0, 1), VECTOR2(0., 1 / sqrt(2)));
}

BOOST_AUTO_TEST_CASE(one_layer) {
  DevicePtr_t device = hpp::pinocchio::unittest::makeDevice(
      hpp::pinocchio::unittest::HumanoidSimple);
  BOOST_REQUIRE(device);
  device->rootJoint()->lowerBound(0, -1);
  device->rootJoint()->lowerBound(1, -1);
  device->rootJoint()->lowerBound(2, -1);
  device->rootJoint()->upperBound(0, 1);
  device->rootJoint()->upperBound(1, 1);
  device->rootJoint()->upperBound(2, 1);
  JointPtr_t ee1 = device->getJointByName("lleg5_joint"),
             ee2 = device->getJointByName("rleg5_joint");

  Configuration_t q = device->currentConfiguration(),
                  qrand = ::pinocchio::randomConfiguration(device->model());

  solver::HierarchicalIterative solver(device->configSpace());
  solver.maxIterations(20);
  solver.errorThreshold(1e-3);
  solver.saturation(hpp::make_shared<solver::saturation::Device>(device));

  device->currentConfiguration(q);
  device->computeForwardKinematics();
  Transform3f tf1(ee1->currentTransformation());
  Transform3f tf2(ee2->currentTransformation());

  ImplicitPtr_t constraint(Implicit::create(
      Orientation::create("Orientation", device, ee2, tf2), 3 * Equality));
  BOOST_CHECK(constraint->comparisonType() == 3 * Equality);
  solver.add(constraint, 0);
  constraint = Implicit::create(Position::create("Position", device, ee2, tf2),
                                3 * Equality);
  BOOST_CHECK(constraint->comparisonType() == 3 * Equality);
  solver.add(constraint, 0);

  BOOST_CHECK(solver.numberStacks() == 1);

  BOOST_CHECK(solver.isSatisfied(q));

  Configuration_t tmp = qrand;
  BOOST_CHECK_EQUAL(solver.solve<solver::lineSearch::Backtracking>(qrand),
                    solver::HierarchicalIterative::SUCCESS);
  qrand = tmp;
  BOOST_CHECK_EQUAL(solver.solve<solver::lineSearch::ErrorNormBased>(qrand),
                    solver::HierarchicalIterative::SUCCESS);
  qrand = tmp;
  BOOST_CHECK_EQUAL(solver.solve<solver::lineSearch::FixedSequence>(qrand),
                    solver::HierarchicalIterative::SUCCESS);
}

template <typename LineSearch = solver::lineSearch::Constant>
struct test_affine_opt : test_base<LineSearch> {
  test_affine_opt(const matrix_t& A, const matrix_t& B)
      : test_base<LineSearch>(A.cols()) {
    // min  X^T * B * X
    // s.t. A * X - 1 = 0
    //      0 <= X <= 1
    BOOST_REQUIRE_EQUAL(A.cols(), B.cols());
    BOOST_REQUIRE_EQUAL(A.rows(), 1);
    BOOST_TEST_MESSAGE(A);
    BOOST_TEST_MESSAGE(B);
    AffineFunctionPtr_t f(AffineFunction::create(A, vector_t::Constant(1, -1)));
    Quadratic::Ptr_t cost(new Quadratic(B));
    ImplicitPtr_t f_constraint(Implicit::create(
        f, ComparisonTypes_t(f->outputDerivativeSize(), Equality)));
    ImplicitPtr_t cost_constraint(Implicit::create(
        cost, ComparisonTypes_t(f->outputDerivativeSize(), Equality)));
    this->solver.add(f_constraint, 0);
    this->solver.add(cost_constraint, 1);
    // this->solver.add(cost, 0);
    this->solver.lastIsOptional(true);
    BOOST_CHECK(this->solver.numberStacks() == 2);
  }

  vector_t optimize(value_type x0, value_type x1) {
    vector_t x(VECTOR2(x0, x1));
    this->solver.lastIsOptional(false);
    this->solver.solve(x, this->ls);
    this->solver.lastIsOptional(true);
    return x;
  }
};

BOOST_AUTO_TEST_CASE(affine_opt) {
  matrix_t A(1, 2);
  A << 1, 1;
  matrix_t B(2, 2);
  B << 1, 0, 0, 1;

  test_affine_opt<> test(A, B);
  test.success(0, 0);
  test.success(0.1, 0);
  test.success(0, 0.1);
  test.success(0.5, 0.5);

  EIGEN_VECTOR_IS_APPROX(test.optimize(0.1, 0), VECTOR2(0.5, 0.5));
  EIGEN_VECTOR_IS_APPROX(test.optimize(0, 0.1), VECTOR2(0.5, 0.5));
  EIGEN_VECTOR_IS_APPROX(test.optimize(0.5, 0.5), VECTOR2(0.5, 0.5));
}

// build an implicit constraint with values in SE3 and with non trivial mask
BOOST_AUTO_TEST_CASE(mask) {
  struct Identity : public DifferentiableFunction {
    static DifferentiableFunctionPtr_t create() {
      return DifferentiableFunctionPtr_t(new Identity());
    }
    Identity() : DifferentiableFunction(7, 6, LiegroupSpace::R3xSO3()) {}
    virtual void impl_compute(LiegroupElementRef result,
                              vectorIn_t argument) const {
      result.vector() = argument;
    }

    virtual void impl_jacobian(matrixOut_t jacobian, vectorIn_t) const {
      jacobian.setIdentity();
    }

    bool isEqual(const DifferentiableFunction& other) const {
      dynamic_cast<const Identity&>(other);
      if (!DifferentiableFunction::isEqual(other)) return false;

      return true;
    }
  };  // class Identity
  solver::HierarchicalIterative solver(LiegroupSpace::R3xSO3());
  solver.maxIterations(20);
  solver.errorThreshold(1e-10);
  std::vector<bool> mask{true, true, false, false, false, true};
  ImplicitPtr_t c1(Implicit::create(Identity::create(), 6 * Equality, mask));
  ImplicitPtr_t c2(Implicit::create(Identity::create(), 6 * Equality, mask));
  solver.add(c1, 0);
  try {
    solver.add(c2, 0);
    BOOST_CHECK(false);
  } catch (const std::logic_error& err) {
    BOOST_CHECK(std::string(err.what()) ==
                std::string("Contraint \"\" already in solver"));
  }
  vector_t q(7);
  q << 1, 2, 3, .5, .5, .5, .5;
  solver.rightHandSideFromConfig(q);
  bool found;
  vector_t error(6);
  std::cout << "q=" << q.transpose() << std::endl;
  BOOST_CHECK(solver.isSatisfied(q));
  BOOST_CHECK(solver.isConstraintSatisfied(c1, q, error, found));
  std::cout << "error=" << error.transpose() << std::endl;
  BOOST_CHECK(found);
  BOOST_CHECK(error.norm() < 1e-10);
  std::cout << solver << std::endl;
  q << 0, 0, 0, 0, 0, 0, 1;
  std::cout << "q=" << q.transpose() << std::endl;
  BOOST_CHECK(!solver.isConstraintSatisfied(c1, q, error, found));
  std::cout << "error=" << error.transpose() << std::endl;
}
