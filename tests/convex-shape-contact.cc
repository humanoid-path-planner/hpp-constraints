// Copyright (c) 2020, LAAS-CNRS
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

#define BOOST_TEST_MODULE hpp_constraints
#include <../tests/convex-shape-contact-function.hh>
#include <../tests/util.hh>
#include <boost/test/included/unit_test.hpp>
#include <hpp/constraints/convex-shape-contact.hh>
#include <hpp/constraints/explicit/convex-shape-contact.hh>
#include <hpp/constraints/solver/by-substitution.hh>
#include <hpp/pinocchio/liegroup-space.hh>
#include <hpp/pinocchio/urdf/util.hh>
#include <pinocchio/algorithm/joint-configuration.hpp>

using hpp::pinocchio::BodyPtr_t;
using hpp::pinocchio::Configuration_t;
using hpp::pinocchio::ConfigurationPtr_t;
using hpp::pinocchio::Device;
using hpp::pinocchio::DevicePtr_t;
using hpp::pinocchio::JointPtr_t;
using hpp::pinocchio::JointVector_t;
using hpp::pinocchio::LiegroupElement;
using hpp::pinocchio::LiegroupSpace;
using hpp::pinocchio::urdf::loadModelFromString;
using namespace hpp::constraints;

LiegroupSpacePtr_t SE3(LiegroupSpace::SE3());
// This test builds a robot with two freeflyer joints j1 and j2
// To each joint, a box of size (2,2,2) is attached: box1 and box2.
// Two floor contact surfaces are attached to box 1
//  - one on face z=1,
//  - one on face z=-1
// Two object contact surfaces are attached to box2
//  - one on face z=1,
//  - one on face z=-1
// 4 configurations q_init, each one corresponding to a type of contact are
// produced.
// A solver of type BySubstitution is built with only one explicit contact
// constraint of type explicit_::ConvexShapeContact
//
// N random configurations are sampled. For each of them, and for each q_init,
// the solver right hand side is initialized with q_init, the solver solves from
// the random configuration.
// the relative position of box2 with respect to box1 is checked to be the same
// as in configuration q_init.
BOOST_AUTO_TEST_CASE(convexShapeContact) {
  const std::string model(
      "<robot name=\"box\">"
      "  <link name=\"baselink\">"
      "    <collision>"
      "      <origin rpy=\"0 0 0\" xyz=\"0 0 0\"/>"
      "      <geometry>"
      "        <box size=\"2 2 2\"/>"
      "      </geometry>"
      "    </collision>"
      "  </link>"
      "</robot>");

  DevicePtr_t robot(Device::create("two-boxes"));
  // Create two freeflying boxes
  loadModelFromString(robot, 0, "1/", "freeflyer", model, "");
  loadModelFromString(robot, 0, "2/", "freeflyer", model, "");
  assert(robot->nbJoints() == 2);
  JointPtr_t j1(robot->jointAt(0));
  JointPtr_t j2(robot->jointAt(1));
  // Set joint bounds
  for (std::size_t i = 0; i < 2; ++i) {
    vector_t l(7);
    l << -2, -2, -2, -1, -1, -1, -1;
    vector_t u(7);
    u << 2, 2, 2, 1, 1, 1, 1;
    robot->jointAt(i)->lowerBounds(l);
    robot->jointAt(i)->upperBounds(u);
  }
  vector3_t v;
  JointAndShape_t surface;
  JointAndShapes_t surfacesObject1, surfacesObject2;
  // Upper face of box1
  surface.first = j1;
  v << -1., -1., 1;
  surface.second.push_back(v);
  v << 1., -1., 1.;
  surface.second.push_back(v);
  v << 1., 1., 1.;
  surface.second.push_back(v);
  v << -1., 1., 1.;
  surface.second.push_back(v);
  surfacesObject1.push_back(surface);
  // Lower face of box1
  surface.second.clear();
  v << -1., 1., -1.;
  surface.second.push_back(v);
  v << 1., 1., -1.;
  surface.second.push_back(v);
  v << 1., -1., -1.;
  surface.second.push_back(v);
  v << -1., -1., -1;
  surface.second.push_back(v);
  surfacesObject1.push_back(surface);
  // Upper face of box2
  surface.second.clear();
  surface.first = j2;
  v << -1., -1., 1;
  surface.second.push_back(v);
  v << 1., -1., 1.;
  surface.second.push_back(v);
  v << 1., 1., 1.;
  surface.second.push_back(v);
  v << -1., 1., 1.;
  surface.second.push_back(v);
  surfacesObject2.push_back(surface);
  // Lower face of box2
  surface.second.clear();
  v << -1., 1., -1.;
  surface.second.push_back(v);
  v << 1., 1., -1.;
  surface.second.push_back(v);
  v << 1., -1., -1.;
  surface.second.push_back(v);
  v << -1., -1., -1;
  surface.second.push_back(v);
  surfacesObject2.push_back(surface);
  ExplicitPtr_t contactConstraint(explicit_::ConvexShapeContact::create(
      "box1/box2", robot, surfacesObject1, surfacesObject2, 0));
  ConvexShapeContactHoldPtr_t f(HPP_DYNAMIC_PTR_CAST(
      ConvexShapeContactHold, contactConstraint->functionPtr()));
  value_type M(f->contactConstraint()->radius());
  BOOST_CHECK(M > sqrt(2));
  // check that the two joints involved can be retrieved correctly
  std::pair<JointConstPtr_t, JointConstPtr_t> joints =
      f->dependsOnRelPoseBetween(nullptr);
  BOOST_CHECK_EQUAL(Joint::index(joints.first), Joint::index(j1));
  BOOST_CHECK_EQUAL(Joint::index(joints.second), Joint::index(j2));

  joints = f->contactConstraint()->dependsOnRelPoseBetween(nullptr);
  BOOST_CHECK_EQUAL(Joint::index(joints.first), Joint::index(j1));
  BOOST_CHECK_EQUAL(Joint::index(joints.second), Joint::index(j2));

  joints = f->complement()->dependsOnRelPoseBetween(nullptr);
  BOOST_CHECK_EQUAL(Joint::index(joints.first), Joint::index(j1));
  BOOST_CHECK_EQUAL(Joint::index(joints.second), Joint::index(j2));
  // box 1 above box 2: surfaces in contact are
  //  - floor surface 1 for box 1,
  //  - object surface 0 for box 2.
  std::vector<Configuration_t> q_init;
  Configuration_t q(14);
  // box 1 above box 2: surfaces in contact are
  //  - floor surface 1 for box 1,
  //  - object surface 0 for box 2.
  q << 0, 0, 0, 0, 0, 0, 1, .5, .8, 2, 0, 0, 0, 1;
  q_init.push_back(q);
  // box 1 above box 2: surfaces in contact are
  //  - floor surface 1 for box 1,
  //  - object surface 1 for box 2.
  q << 0, 0, 0, 0, 0, 0, 1, -.25, .2, 2, 1, 0, 0, 0;
  q_init.push_back(q);
  // box 1 below box 2: surfaces in contact are
  //  - floor surface 0 for box 1,
  //  - object surface 1 for box 2.
  q << 0, 0, 0, 0, 0, 0, 1, .2, .35, -2, 0, 0, 0, 1;
  q_init.push_back(q);
  // box 1 below box 2: surfaces in contact are
  //  - floor surface 0 for box 1,
  //  - object surface 0 for box 2.
  q << 0, 0, 0, 0, 0, 0, 1, -.8, .6, -2, 1, 0, 0, 0;
  q_init.push_back(q);
  std::vector<LiegroupElement> q1Invq2Exp;
  // Create Solver by substitution with one explicit constraint.
  const value_type epsilon(1e-7);
  solver::BySubstitution solver1(SE3 * SE3);
  solver::HierarchicalIterative solver2(SE3 * SE3);
  solver1.errorThreshold(epsilon);
  solver2.errorThreshold(epsilon);
  solver2.maxIterations(30);
  solver1.add(contactConstraint);
  solver2.add(contactConstraint, 0);
  vector_t error(8);
  // Compute expected box relative positions for each q_init
  for (std::size_t i = 0; i < 4; ++i) {
    LiegroupElement q1(q_init[i].head<7>(), SE3);
    LiegroupElement q2(q_init[i].tail<7>(), SE3);
    // log of position of box2 with respect to box1.
    q1Invq2Exp.push_back(SE3->exp(q2 - q1));
    // Check that q_init satisfies the constraint
    vector_t rhs(solver1.rightHandSideFromConfig(q_init[i]));
    BOOST_CHECK(solver1.isSatisfied(q_init[i], error));
    rhs = solver2.rightHandSideFromConfig(q_init[i]);
    BOOST_CHECK(solver2.isSatisfied(q_init[i]));
  }
  const std::size_t N(100);
  // Generate random configurations
  std::vector<Configuration_t> q_rand;
  for (std::size_t i = 0; i < N; ++i) {
    q_rand.push_back(::pinocchio::randomConfiguration(robot->model()));
  }
  q_rand[0].head<6>().setZero();
  q_rand[0][6] = 1;

  std::size_t nSuccesses(0);
  for (std::size_t i = 0; i < N; ++i) {
    Configuration_t q(q_rand[i]);
    for (std::size_t j = 0; j < 4; ++j) {
      solver1.rightHandSideFromConfig(q_init[j]);
      solver2.rightHandSideFromConfig(q_init[j]);
      BOOST_CHECK(solver2.isSatisfied(q_init[j]));
      q = q_rand[i];
      solver::HierarchicalIterative::Status status(solver1.solve(q));
      BOOST_CHECK(status == solver::HierarchicalIterative::SUCCESS);
      // Check that box 1 did not move.
      BOOST_CHECK(q.head<7>() == q_rand[i].head<7>());
      LiegroupElement q1(q.head<7>(), SE3);
      LiegroupElement q2(q.tail<7>(), SE3);
      LiegroupElement q1Invq2(SE3->exp(q2 - q1));
      // Check that relative position of box 2 with respect to box 1
      // is the same as defined by q_init.
      BOOST_CHECK((q1Invq2 - q1Invq2Exp[j]).norm() <=
                  10 * solver1.errorThreshold());
      q = q_rand[i];
      status = solver2.solve(q);
      if (status == solver::HierarchicalIterative::SUCCESS) {
        ++nSuccesses;
        q1.vector() = q.head<7>();
        q2.vector() = q.tail<7>();
        q1Invq2 = SE3->exp(q2 - q1);
        // Note that error of constraint is not equivalent to relative
        // position between objects, but between contact surface frames.
        BOOST_CHECK((q1Invq2 - q1Invq2Exp[j]).norm() <=
                    10 * solver2.errorThreshold());
      } else {
        solver2.isSatisfied(q);
        solver2.residualError(error);
      }
    }
  }
  // Check that success rate is not too low. N/10 is an arbitrary value.
  BOOST_CHECK(nSuccesses >= N / 10);
}
