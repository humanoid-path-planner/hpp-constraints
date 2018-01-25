// Copyright (c) 2017, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//
// This file is part of hpp-constraints.
// hpp-constraints is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-constraints is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-constraints. If not, see <http://www.gnu.org/licenses/>.

#define BOOST_TEST_MODULE HIERARCHICAL_ITERATIVE_SOLVER
#include <boost/test/unit_test.hpp>

#include <hpp/constraints/iterative-solver.hh>

#include <functional>

#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/simple-device.hh>
#include <hpp/constraints/generic-transformation.hh>

#include <../tests/util.hh>

using namespace hpp::constraints;
value_type test_precision = 1e-6;

#define EIGEN_IS_APPROX(matrixA, matrixB)                                      \
  BOOST_CHECK_MESSAGE(matrixA.isApprox(matrixB, test_precision),               \
      "check " #matrixA ".isApprox(" #matrixB ") failed "                      \
      "[" << matrixA.transpose() << " != " << matrixB.transpose() << "]")

#define EIGEN_IS_NOT_APPROX(matrixA, matrixB)                                  \
  BOOST_CHECK_MESSAGE(!matrixA.isApprox(matrixB, test_precision),              \
      "check !" #matrixA ".isApprox(" #matrixB ") failed "                     \
      "[" << matrixA.transpose() << " != " << matrixB.transpose() << "]")

#define VECTOR2(x0, x1) ((hpp::constraints::vector_t (2) << x0, x1).finished())

using Eigen::VectorXi;

void addition (vectorIn_t from, vectorIn_t velocity, vectorOut_t result)
{
  // result = from + velocity;
  result = (from + velocity).cwiseMin(1).cwiseMax(-1);
}
bool saturation (vectorIn_t x, VectorXi& sat)
{
  bool ret = false;
  for (size_type i = 0; i < x.size(); ++i) {
    if (x[i] <= 0) {
      sat[i] = -1;
      ret = true;
    }
    else if (x[i] >= 1) {
      sat[i] =  1;
      ret = true;
    } else {
      sat[i] = 0;
    }
  }
  return ret;
}

template <typename LineSearch = lineSearch::Constant>
struct test_quadratic
{
  HierarchicalIterativeSolver solver;
  LineSearch ls;

  test_quadratic (const matrix_t& A) : solver(2, 2)
  {
    // Find (x, y)
    // s.t. a * x^2 + b * y^2 - 1 = 0
    //      0 <= x <= 1
    //      0 <= y <= 1
    BOOST_TEST_MESSAGE(A);
    Quadratic::Ptr_t f (new Quadratic (A, -1));

    solver.maxIterations(20);
    solver.errorThreshold(test_precision);
    solver.integration(addition);
    solver.saturation(saturation);

    solver.add(f, 0);
    BOOST_CHECK(solver.numberStacks() == 1);
  }

  vector_t success (value_type x0, value_type x1)
  {
    vector_t x (VECTOR2(x0,x1));
    BOOST_CHECK_EQUAL(solver.solve(x, ls), HierarchicalIterativeSolver::SUCCESS);
    return x;
  }

  vector_t failure (value_type x0, value_type x1)
  {
    vector_t x (VECTOR2(x0,x1));
    BOOST_CHECK_PREDICATE (std::not_equal_to<HierarchicalIterativeSolver::Status>(), (solver.solve(x, ls))(HierarchicalIterativeSolver::SUCCESS));
    return x;
  }
};

BOOST_AUTO_TEST_CASE(quadratic)
{
  matrix_t A(2,2);

  A << 1, 0, 0, 1;
  test_quadratic<> test (A);
  BOOST_CHECK_EQUAL (test.failure(0,0), VECTOR2(0,0));
  test.success(0.1,0);
  test.success(0,0.1);
  test.success(0.5, 0.5);

  A << 2, 0, 0, 2;
  test = test_quadratic<> (A);
  test.success(0.1,0);
  test.success(0,0.1);
  test.success(0.5, 0.5);

  A << 0.5, 0, 0, 0.5;
  test = test_quadratic<> (A);
  // This is exact because of the saturation
  BOOST_CHECK_EQUAL (test.success (1, 0.001), VECTOR2(1,1)); // Slide on the border x = 1
  BOOST_CHECK_EQUAL (test.success (0.001, 1), VECTOR2(1,1)); // Slide on the border y = 1

  A << 0.75, 0, 0, 0.75;
  test_quadratic<lineSearch::FixedSequence> test4 (A);
  // This is not exact because the solver does not saturate.
  EIGEN_IS_APPROX (test4.success (1, 0.1), VECTOR2(1.,1/sqrt(3))); // Slide on the border x = 1
  EIGEN_IS_APPROX (test4.success (0.1, 1), VECTOR2(1/sqrt(3),1.)); // Slide on the border y = 1
  // There is an overshoot. To overcome this, the Hessian of the function should be obtained.
  EIGEN_IS_NOT_APPROX (test4.success (1, 0.001), VECTOR2(1.,1/sqrt(3))); // Slide on the border x = 1
  EIGEN_IS_NOT_APPROX (test4.success (0.001, 1), VECTOR2(1/sqrt(3),1.)); // Slide on the border y = 1

  // Ellipsoid: computations are approximative
  A << 0.5, 0, 0, 2;
  test_quadratic<lineSearch::FixedSequence> test1 (A);
  BOOST_CHECK_EQUAL (test1.success (1, 0.5), VECTOR2(1.,0.5)); // Slide on the border x = 1
  EIGEN_IS_APPROX (test1.success (1, 0.1), VECTOR2(1.,0.5)); // Slide on the border x = 1
  EIGEN_IS_APPROX (test1.success (0, 1), VECTOR2(0.,1/sqrt(2)));
}

BOOST_AUTO_TEST_CASE(one_layer)
{
  DevicePtr_t device = hpp::pinocchio::unittest::makeDevice (hpp::pinocchio::unittest::HumanoidRomeo);
  BOOST_REQUIRE (device);
  device->rootJoint()->lowerBound (0, -1);
  device->rootJoint()->lowerBound (1, -1);
  device->rootJoint()->lowerBound (2, -1);
  device->rootJoint()->upperBound (0,  1);
  device->rootJoint()->upperBound (1,  1);
  device->rootJoint()->upperBound (2,  1);
  JointPtr_t ee1 = device->getJointByName ("LAnkleRoll"),
             ee2 = device->getJointByName ("RAnkleRoll");

  Configuration_t q = device->currentConfiguration (),
                  qrand = se3::randomConfiguration(device->model());

  HierarchicalIterativeSolver solver(device->configSize(), device->numberDof());
  solver.maxIterations(20);
  solver.errorThreshold(1e-3);
  solver.integration(boost::bind(hpp::pinocchio::integrate<true, se3::LieGroupTpl>, device, _1, _2, _3));
  solver.saturation(boost::bind(saturate, device, _1, _2));

  device->currentConfiguration (q);
  device->computeForwardKinematics ();
  Transform3f tf1 (ee1->currentTransformation ());
  Transform3f tf2 (ee2->currentTransformation ());

  solver.add(Orientation::create ("Orientation", device, ee2, tf2), 0);
  solver.add(Position::create    ("Position"   , device, ee2, tf2), 0);

  BOOST_CHECK(solver.numberStacks() == 1);

  BOOST_CHECK(solver.isSatisfied(q));

  Configuration_t tmp = qrand;
  BOOST_CHECK_EQUAL(solver.solve<lineSearch::Backtracking  >(qrand), HierarchicalIterativeSolver::SUCCESS);
  qrand = tmp;
  BOOST_CHECK_EQUAL(solver.solve<lineSearch::ErrorNormBased>(qrand), HierarchicalIterativeSolver::SUCCESS);
  qrand = tmp;
  BOOST_CHECK_EQUAL(solver.solve<lineSearch::FixedSequence >(qrand), HierarchicalIterativeSolver::SUCCESS);
}

