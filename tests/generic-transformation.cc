// Copyright (c) 2016, Joseph Mirabel
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

#define EIGEN_RUNTIME_NO_MALLOC

#include <hpp/constraints/explicit/relative-pose.hh>
#include <hpp/constraints/generic-transformation.hh>
#include <hpp/constraints/solver/by-substitution.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint-collection.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/serialization.hh>
#include <hpp/pinocchio/simple-device.hh>
#include <hpp/pinocchio/urdf/util.hh>
#include <pinocchio/algorithm/joint-configuration.hpp>
#include <sstream>

#include "hpp/constraints/tools.hh"

#define BOOST_TEST_MODULE hpp_constraints
#include <stdlib.h>

#include <boost/test/included/unit_test.hpp>

using hpp::constraints::EqualToZero;
using hpp::pinocchio::Configuration_t;
using hpp::pinocchio::ConfigurationPtr_t;
using hpp::pinocchio::Device;
using hpp::pinocchio::DevicePtr_t;
using hpp::pinocchio::JOINT_POSITION;
using hpp::pinocchio::JointIndex;
using hpp::pinocchio::JointPtr_t;
using hpp::pinocchio::LiegroupSpace;
using hpp::pinocchio::Transform3f;

using hpp::pinocchio::urdf::loadModelFromString;

using namespace hpp::constraints;

class BasicConfigurationShooter {
 public:
  BasicConfigurationShooter(const DevicePtr_t& robot) : robot_(robot) {}
  virtual Configuration_t shoot() const {
    size_type extraDim = robot_->extraConfigSpace().dimension();
    size_type offset = robot_->configSize() - extraDim;

    Configuration_t config(robot_->configSize());
    config.head(offset) = ::pinocchio::randomConfiguration(robot_->model());

    // Shoot extra configuration variables
    for (size_type i = 0; i < extraDim; ++i) {
      value_type lower = robot_->extraConfigSpace().lower(i);
      value_type upper = robot_->extraConfigSpace().upper(i);
      value_type range = upper - lower;
      if ((range < 0) || (range == std::numeric_limits<double>::infinity())) {
        std::ostringstream oss("Cannot uniformy sample extra config variable ");
        oss << i << ". min = " << lower << ", max = " << upper << std::endl;
        throw std::runtime_error(oss.str());
      }
      config[offset + i] = lower + (upper - lower) * rand() / RAND_MAX;
    }
    return config;
  }

 private:
  DevicePtr_t robot_;
};  // class BasicConfigurationShooter

BOOST_AUTO_TEST_CASE(print) {
  DevicePtr_t device = hpp::pinocchio::unittest::makeDevice(
      hpp::pinocchio::unittest::HumanoidSimple);
  JointPtr_t ee1 = device->getJointByName("lleg5_joint"),
             ee2 = device->getJointByName("rleg5_joint");
  BOOST_REQUIRE(device);
  BasicConfigurationShooter cs(device);

  device->currentConfiguration(cs.shoot());
  device->computeForwardKinematics(JOINT_POSITION);
  Transform3f tf1(ee1->currentTransformation());
  Transform3f tf2(ee2->currentTransformation());

  std::vector<DifferentiableFunctionPtr_t> functions;
  functions.push_back(Orientation::create("Orientation", device, ee2, tf2));
  functions.push_back(Position::create("Position", device, ee2, tf2, tf1));
  functions.push_back(
      Transformation::create("Transformation", device, ee1, tf1));
  functions.push_back(RelativeOrientation::create("RelativeOrientation", device,
                                                  ee1, ee2, tf1));
  functions.push_back(
      RelativePosition::create("RelativePosition", device, ee1, ee2, tf1, tf2));
  functions.push_back(RelativeTransformation::create(
      "RelativeTransformation", device, ee1, ee2, tf1, tf2));

  Configuration_t q1 = cs.shoot(), q2 = cs.shoot();
  for (std::size_t i = 0; i < functions.size(); ++i) {
    DifferentiableFunctionPtr_t f = functions[i];

    std::cout << *f << std::endl;

    LiegroupElement v(f->outputSpace());
    matrix_t J(f->outputDerivativeSize(), f->inputDerivativeSize());

    f->value(v, q1);
    f->jacobian(J, q1);
    // TODO this is broken at the moment because of the introduction
    // of a multithreaded device.
    // Eigen::internal::set_is_malloc_allowed(false);
    f->value(v, q2);
    f->jacobian(J, q2);
    // Eigen::internal::set_is_malloc_allowed(true);
  }

  // Check active parameters
  ArrayXb ap1 =
      Orientation::create("Orientation", device, ee1, tf1)->activeParameters();
  ArrayXb ap2 =
      Orientation::create("Orientation", device, ee2, tf2)->activeParameters();
  ArrayXb ap12 =
      RelativeOrientation::create("RelativeOrientation", device, ee1, ee2, tf1)
          ->activeParameters();

  ArrayXb not_ap1 = (ap1 == false);
  ArrayXb not_ap2 = (ap2 == false);
  BOOST_CHECK((ap12 == ((not_ap1 && ap2) || (ap1 && not_ap2))).all());

  // Check active derivative parameters
  ap1 = Orientation::create("Orientation", device, ee1, tf1)
            ->activeDerivativeParameters();
  ap2 = Orientation::create("Orientation", device, ee2, tf2)
            ->activeDerivativeParameters();
  ap12 =
      RelativeOrientation::create("RelativeOrientation", device, ee1, ee2, tf1)
          ->activeDerivativeParameters();

  not_ap1 = (ap1 == false);
  not_ap2 = (ap2 == false);
  BOOST_CHECK((ap12 == ((not_ap1 && ap2) || (ap1 && not_ap2))).all());
}

BOOST_AUTO_TEST_CASE(multithread) {
  DevicePtr_t device = hpp::pinocchio::unittest::makeDevice(
      hpp::pinocchio::unittest::HumanoidSimple);
  device->numberDeviceData(4);
  JointPtr_t ee1 = device->getJointByName("lleg5_joint"),
             ee2 = device->getJointByName("rleg5_joint");
  BOOST_REQUIRE(device);
  BasicConfigurationShooter cs(device);

  device->currentConfiguration(cs.shoot());
  device->computeForwardKinematics(JOINT_POSITION);
  Transform3f tf1(ee1->currentTransformation());
  Transform3f tf2(ee2->currentTransformation());

  std::vector<DifferentiableFunctionPtr_t> functions;
  functions.push_back(Orientation::create("Orientation", device, ee2, tf2));
  functions.push_back(Position::create("Position", device, ee2, tf2, tf1));
  functions.push_back(
      Transformation::create("Transformation", device, ee1, tf1));
  functions.push_back(RelativeOrientation::create("RelativeOrientation", device,
                                                  ee1, ee2, tf1));
  functions.push_back(
      RelativePosition::create("RelativePosition", device, ee1, ee2, tf1, tf2));
  functions.push_back(RelativeTransformation::create(
      "RelativeTransformation", device, ee1, ee2, tf1, tf2));
  functions.push_back(RelativeOrientation::create("RelativeOrientation", device,
                                                  ee1, JointPtr_t(), tf1));
  functions.push_back(RelativePosition::create("RelativePosition", device, ee1,
                                               JointPtr_t(), tf1, tf2));
  functions.push_back(RelativeTransformation::create(
      "RelativeTransformation", device, ee1, JointPtr_t(), tf1, tf2));

  const int N = 10;
  Configuration_t q = cs.shoot();
  for (std::size_t i = 0; i < functions.size(); ++i) {
    DifferentiableFunctionPtr_t f = functions[i];

    std::vector<LiegroupElement> vs(N, LiegroupElement(f->outputSpace()));
    std::vector<matrix_t> Js(
        N, matrix_t(f->outputDerivativeSize(), f->inputDerivativeSize()));
#pragma omp parallel for
    for (int j = 0; j < 10; ++j) {
      f->value(vs[j], q);
      f->jacobian(Js[j], q);
    }

    for (int j = 1; j < N; ++j) {
      BOOST_CHECK_EQUAL(vs[0].vector(), vs[j].vector());
      BOOST_CHECK_EQUAL(Js[0], Js[j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(serialization) {
  DevicePtr_t device = hpp::pinocchio::unittest::makeDevice(
      hpp::pinocchio::unittest::HumanoidSimple);
  device->numberDeviceData(4);
  JointPtr_t ee1 = device->getJointByName("lleg5_joint"),
             ee2 = device->getJointByName("rleg5_joint");
  BOOST_REQUIRE(device);

  device->currentConfiguration(device->neutralConfiguration());
  device->computeForwardKinematics(JOINT_POSITION);
  Transform3f tf1(ee1->currentTransformation());
  Transform3f tf2(ee2->currentTransformation());

  std::vector<DifferentiableFunctionPtr_t> functions;
  functions.push_back(Orientation::create("Orientation", device, ee2, tf2));
  functions.push_back(Position::create("Position", device, ee2, tf2, tf1));
  functions.push_back(
      Transformation::create("Transformation", device, ee1, tf1));
  functions.push_back(RelativeOrientation::create("RelativeOrientation", device,
                                                  ee1, ee2, tf1));
  functions.push_back(
      RelativePosition::create("RelativePosition", device, ee1, ee2, tf1, tf2));
  functions.push_back(RelativeTransformation::create(
      "RelativeTransformation", device, ee1, ee2, tf1, tf2));
  functions.push_back(RelativeOrientation::create("RelativeOrientation", device,
                                                  ee1, JointPtr_t(), tf1));
  functions.push_back(RelativePosition::create("RelativePosition", device, ee1,
                                               JointPtr_t(), tf1, tf2));
  functions.push_back(RelativeTransformation::create(
      "RelativeTransformation", device, ee1, JointPtr_t(), tf1, tf2));

  for (std::size_t i = 0; i < functions.size(); ++i) {
    DifferentiableFunctionPtr_t f = functions[i];

    DifferentiableFunctionPtr_t f_restored;

    std::stringstream ss;
    {
      hpp::serialization::xml_oarchive oa(ss);
      oa.insert(device->name(), device.get());
      oa << boost::serialization::make_nvp("function", f);
    }

    {
      hpp::serialization::xml_iarchive ia(ss);
      ia.insert(device->name(), device.get());
      ia >> boost::serialization::make_nvp("function", f_restored);
    }

    std::ostringstream oss1, oss2;
    oss1 << *f;
    oss2 << *f_restored;

    BOOST_CHECK_EQUAL(oss1.str(), oss2.str());
  }
}

BOOST_AUTO_TEST_CASE(RelativeTransformation_R3xSO3) {
  const std::string model(
      "<robot name=\"box\">"
      "  <link name=\"baselink\">"
      "  </link>"
      "</robot>");

  DevicePtr_t robot(Device::create("two-freeflyers"));
  // Create two freeflying boxes
  loadModelFromString(robot, 0, "1/", "freeflyer", model, "");
  loadModelFromString(robot, 0, "2/", "freeflyer", model, "");
  BOOST_CHECK(robot->configSize() == 14);
  BOOST_CHECK(robot->numberDof() == 12);
  BOOST_CHECK(robot->nbJoints() == 2);
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
  // Create constraint
  //
  // Excerpt from romeo-placard benchmark
  // Joint1: romeo/LWristPitch
  // Frame in joint 1
  //   R = 0.7071067739978436073, 0.70710678837525142715, 0
  //       -2.2663502965461253728e-09,  2.2663502504650490188e-09, -1
  //       -0.70710678837525142715,
  //       0.70710677399784382935, 3.2051032938795742666e-09
  //   p = 0.099999999776482578762, -3.2051032222399330291e-11,
  //   -0.029999999776482582509
  // Joint2: placard/root_joint
  // Frame in joint 2
  //   R =   1, 0, 0
  //         0, 1, 0
  //         0, 0, 1
  //   p = 0, 0, -0.34999999403953552246
  // mask: 1, 1, 1, 1, 1, 1,
  // Rhs: 0, 0, 0, 0, 0, 0, 1
  // active rows: [ 0, 5],

  matrix3_t R1, R2;
  vector3_t p1, p2;
  R1 << 0.7071067739978436073, 0.70710678837525142715, 0,
      -2.2663502965461253728e-09, 2.2663502504650490188e-09, -1,
      -0.70710678837525142715, 0.70710677399784382935,
      3.2051032938795742666e-09;
  p1 << 0.099999999776482578762, -3.2051032222399330291e-11,
      -0.029999999776482582509;
  R2.setIdentity();
  p2 << 0, 0, -0.34999999403953552246;
  Transform3f tf1(R1, p1), tf2(R2, p2);
  std::vector<bool> mask = {false, false, false, false, false, true};
  ImplicitPtr_t constraint(Implicit::create(
      RelativeTransformationR3xSO3::create("RelativeTransformationR3xSO3",
                                           robot, j1, j2, tf1, tf2),
      6 * Equality, mask));
  BasicConfigurationShooter cs(robot);
  solver::BySubstitution solver(robot->configSpace());
  solver.errorThreshold(1e-10);
  solver.add(constraint);
  // Check that after setting right hand side with a configuration
  // the configuration satisfies the constraint since comparison type is
  // Equality.
  for (size_type i = 0; i < 1000; ++i) {
    Configuration_t q(cs.shoot());
    vector6_t error;
    solver.rightHandSideFromConfig(q);
    BOOST_CHECK(solver.isSatisfied(q, error));
  }

  // Create grasp constraint with one degree of freedom in rotation along z
  mask = {true, true, true, true, true, false};
  ImplicitPtr_t c1(Implicit::create(
      RelativeTransformationR3xSO3::create("RelativeTransformationR3xSO3",
                                           robot, j1, j2, tf1, tf2),
      6 * EqualToZero, mask));
  solver::BySubstitution s1(robot->configSpace());
  s1.errorThreshold(1e-10);
  s1.add(c1);
  // Create grasp + complement as an explicit constraint
  ExplicitPtr_t c2(
      explicit_::RelativePose::create("ExplicitRelativePose", robot, j1, j2,
                                      tf1, tf2, 5 * EqualToZero << Equality));
  solver::BySubstitution s2(robot->configSpace());
  s2.errorThreshold(1e-4);
  s2.add(c2);

  for (size_type i = 0; i < 0; ++i) {
    Configuration_t q_near(cs.shoot());
    Configuration_t q_new(cs.shoot());
    if (i == 0) {
      // These configuration reproduce a numerical issue encountered with
      // benhmark romeo-placard.
      // If computation was exact, any configuration satisfying c2 should
      // satisfy c1.
      // Configuration q_new below satisfies c2 but not c1.
      q_near << 0.18006349590534418, 0.3627623741970175, 0.9567759630330663,
          0.044416054309488175, 0.31532356328825556, 0.4604329042168087,
          0.8286131819306576, 0.45813483973344404, 0.23514459283216355,
          0.7573015903787429, 0.8141495491430896, 0.1383820163733335,
          0.3806970356973106, 0.4160296818567576;
      q_new << 0.16026892741853033, 0.33925098736439646, 0.8976880203169203,
          -0.040130835169737825, 0.37473431876437147, 0.4405275981290593,
          0.8148000624051422, 0.43787674119234027, 0.18316291571416676,
          0.7189377922181226, 0.7699579340925136, 0.1989432638510445,
          0.35960786236482944, 0.4881275886709128;
    }
    s2.rightHandSideFromConfig(q_near);
    vector6_t error;
    BOOST_CHECK(s1.isSatisfied(q_near, error));
    hppDout(info, error);
    BOOST_CHECK(s2.isSatisfied(q_near, error));
    hppDout(info, error);
    BOOST_CHECK(s1.isSatisfied(q_new, error));
    hppDout(info, error);
    BOOST_CHECK(s2.isSatisfied(q_new, error));
    hppDout(info, error);

    hppDout(info, s1);
    hppDout(info, s2);
  }
}

BOOST_AUTO_TEST_CASE(equality) {
  DevicePtr_t device = hpp::pinocchio::unittest::makeDevice(
      hpp::pinocchio::unittest::HumanoidSimple);
  JointPtr_t ee1 = device->getJointByName("lleg5_joint"),
             ee2 = device->getJointByName("rleg5_joint");
  BOOST_REQUIRE(device);
  BasicConfigurationShooter cs(device);

  device->currentConfiguration(cs.shoot());
  device->computeForwardKinematics(JOINT_POSITION);
  Transform3f tf1(ee1->currentTransformation());
  Transform3f tf2(ee2->currentTransformation());

  std::vector<DifferentiableFunctionPtr_t> functions;
  functions.push_back(Orientation::create("Orientation", device, ee2, tf2));
  functions.push_back(RelativeTransformation::create("Orientation", device, ee1,
                                                     ee2, tf1, tf2));
  functions.push_back(RelativeTransformation::create(
      "RelativeTransformation", device, ee1, ee2, tf1, tf2));
  functions.push_back(RelativeTransformation::create(
      "RelativeTransformation", device, ee1, ee2, tf1, tf2));
  functions.push_back(RelativeTransformation::create(
      "othername_____________", device, ee1, ee2, tf1, tf2));
  // functions[2] and functions[3] are meant to have the same value with
  // different pointers

  BOOST_CHECK(functions[2].get() !=
              functions[3].get());  // boost implementation for ==
  BOOST_CHECK(*functions[2] ==
              *functions[3]);  // uses operator== defined in DiffFunc
  BOOST_CHECK(*functions[0] != *functions[2]);  // a lot of things are different
  BOOST_CHECK(*functions[2] != *functions[4]);  // only the names are different
  BOOST_CHECK(*functions[0] != *functions[1]);  // only the names are equal
}

/* Create a robot with the following kinematic chain. All joints are
   translations along x.

                               universe
                                  |Px
                               test_x
                             /Px       \Px
                       joint_a0       joint_b0
                           |Px            |Px
                       joint_a1       joint_b1
                                          |FF
                                      joint_b2

*/
DevicePtr_t createRobot() {
  std::string urdf(
      "<robot name='test'>"
      "<link name='base_link'/>"
      "<link name='link_test_x'/>"
      "<joint name='test_x' type='prismatic'>"
      "<parent link='base_link'/>"
      "<child  link='link_test_x'/>"
      "<limit effort='30' velocity='1.0' lower='-4' upper='4'/>"
      "</joint>"

      "<link name='link_a0'/>"
      "<link name='link_a1'/>"
      "<joint name='joint_a0' type='prismatic'>"
      "<parent link='link_test_x'/>"
      "<child  link='link_a0'/>"
      "<limit effort='30' velocity='1.0' lower='-4' upper='4'/>"
      "</joint>"
      "<joint name='joint_a1' type='prismatic'>"
      "<parent link='link_a0'/>"
      "<child  link='link_a1'/>"
      "<limit effort='30' velocity='1.0' lower='-4' upper='4'/>"
      "</joint>"

      "<link name='link_b0'/>"
      "<link name='link_b1'/>"
      "<link name='link_b2'/>"
      "<joint name='joint_b0' type='prismatic'>"
      "<parent link='link_test_x'/>"
      "<child  link='link_b0'/>"
      "<limit effort='30' velocity='1.0' lower='-4' upper='4'/>"
      "</joint>"
      "<joint name='joint_b1' type='prismatic'>"
      "<parent link='link_b0'/>"
      "<child  link='link_b1'/>"
      "<limit effort='30' velocity='1.0' lower='-4' upper='4'/>"
      "</joint>"
      "<joint name='joint_b2' type='floating'>"
      "<parent link='link_b1'/>"
      "<child  link='link_b2'/>"
      "</joint>"

      "</robot>");

  DevicePtr_t robot = Device::create("test");
  loadModelFromString(robot, 0, "", "anchor", urdf, "");
  return robot;
}

BOOST_AUTO_TEST_CASE(dependsOnRelPoseBetween) {
  DevicePtr_t device = createRobot();
  BOOST_REQUIRE(device);

  JointPtr_t ee1 = device->getJointByName("joint_a1"),
             ee2 = device->getJointByName("joint_b2");
  BOOST_REQUIRE(device);
  // ensure that the joint indices are as expected
  BOOST_REQUIRE(ee1->index() < ee2->index());

  device->currentConfiguration(device->neutralConfiguration());
  device->computeForwardKinematics(JOINT_POSITION);
  Transform3f tf1(ee1->currentTransformation());
  Transform3f tf2(ee2->currentTransformation());

  DifferentiableFunctionPtr_t function;
  std::pair<JointConstPtr_t, JointConstPtr_t> joints;
  std::pair<JointConstPtr_t, JointConstPtr_t> jointsConstrained;
  ImplicitPtr_t constraint;

  function = Orientation::create("Orientation", device, ee2, tf2);
  joints = function->dependsOnRelPoseBetween(nullptr);
  BOOST_CHECK(!joints.first);
  BOOST_CHECK_EQUAL(joints.second->index(), ee2->index());
  // constraint does not fully constrain the relative pose
  // since it only involves the orientation
  constraint = Implicit::create(
      function, ComparisonTypes_t(function->outputSpace()->nv(), EqualToZero),
      std::vector<bool>(function->outputSpace()->nv(), true));
  jointsConstrained = constraint->doesConstrainRelPoseBetween(device);
  BOOST_CHECK(!jointsConstrained.first);
  BOOST_CHECK(!jointsConstrained.second);

  function = Position::create("Position", device, ee2, tf2, tf1);
  joints = function->dependsOnRelPoseBetween(nullptr);
  BOOST_CHECK(!joints.first);
  BOOST_CHECK_EQUAL(joints.second->index(), ee2->index());
  // constraint does not fully constrain the relative pose
  // since it only involves the position
  constraint = Implicit::create(
      function, ComparisonTypes_t(function->outputSpace()->nv(), EqualToZero),
      std::vector<bool>(function->outputSpace()->nv(), true));
  jointsConstrained = constraint->doesConstrainRelPoseBetween(device);
  BOOST_CHECK(!jointsConstrained.first);
  BOOST_CHECK(!jointsConstrained.second);

  function = Transformation::create("Transformation", device, ee1, tf1);
  joints = function->dependsOnRelPoseBetween(nullptr);
  BOOST_CHECK(!joints.first);
  BOOST_CHECK_EQUAL(joints.second->index(), ee1->index());
  // constraint does not fully constrain the relative pose
  // since the mask is not full
  constraint = Implicit::create(
      function, ComparisonTypes_t(function->outputSpace()->nv(), EqualToZero),
      std::vector<bool>(function->outputSpace()->nv(), false));
  jointsConstrained = constraint->doesConstrainRelPoseBetween(device);
  BOOST_CHECK(!jointsConstrained.first);
  BOOST_CHECK(!jointsConstrained.second);
  // constraint fully constrains the relative pose
  constraint = Implicit::create(
      function, ComparisonTypes_t(function->outputSpace()->nv(), EqualToZero),
      std::vector<bool>(function->outputSpace()->nv(), true));
  jointsConstrained = constraint->doesConstrainRelPoseBetween(device);
  BOOST_CHECK(!jointsConstrained.first);
  BOOST_CHECK_EQUAL(jointsConstrained.second->index(), ee1->index());

  function =
      RelativeOrientation::create("RelativeOrientation", device, ee1, ee2, tf1);
  joints = function->dependsOnRelPoseBetween(nullptr);
  BOOST_CHECK_EQUAL(joints.first->index(), ee1->index());
  BOOST_CHECK_EQUAL(joints.second->index(), ee2->index());
  // constraint does not fully constrain the relative pose
  constraint = Implicit::create(
      function, ComparisonTypes_t(function->outputSpace()->nv(), EqualToZero),
      std::vector<bool>(function->outputSpace()->nv(), true));
  jointsConstrained = constraint->doesConstrainRelPoseBetween(device);
  BOOST_CHECK(!jointsConstrained.first);
  BOOST_CHECK(!jointsConstrained.second);

  function =
      RelativePosition::create("RelativePosition", device, ee1, ee2, tf1, tf2);
  joints = function->dependsOnRelPoseBetween(nullptr);
  BOOST_CHECK_EQUAL(joints.first->index(), ee1->index());
  BOOST_CHECK_EQUAL(joints.second->index(), ee2->index());
  // constraint does not fully constrain the relative pose
  constraint = Implicit::create(
      function, ComparisonTypes_t(function->outputSpace()->nv(), EqualToZero),
      std::vector<bool>(function->outputSpace()->nv(), true));
  jointsConstrained = constraint->doesConstrainRelPoseBetween(device);
  BOOST_CHECK(!jointsConstrained.first);
  BOOST_CHECK(!jointsConstrained.second);

  function = RelativeTransformation::create("RelativeTransformation", device,
                                            ee1, ee2, tf1, tf2);
  joints = function->dependsOnRelPoseBetween(nullptr);
  BOOST_CHECK_EQUAL(joints.first->index(), ee1->index());
  BOOST_CHECK_EQUAL(joints.second->index(), ee2->index());
  // constraint does not fully constrain the relative pose
  // since mask is not full
  constraint = Implicit::create(
      function, ComparisonTypes_t(function->outputSpace()->nv(), EqualToZero),
      std::vector<bool>(function->outputSpace()->nv(), false));
  jointsConstrained = constraint->doesConstrainRelPoseBetween(device);
  BOOST_CHECK(!jointsConstrained.first);
  BOOST_CHECK(!jointsConstrained.second);
  // constraint fully constrains the relative pose
  constraint = Implicit::create(
      function, ComparisonTypes_t(function->outputSpace()->nv(), EqualToZero),
      std::vector<bool>(function->outputSpace()->nv(), true));
  jointsConstrained = constraint->doesConstrainRelPoseBetween(device);
  BOOST_CHECK_EQUAL(jointsConstrained.first->index(), ee1->index());
  BOOST_CHECK_EQUAL(jointsConstrained.second->index(), ee2->index());

  function = RelativeOrientation::create("RelativeOrientation", device, ee1,
                                         JointPtr_t(), tf1);
  joints = function->dependsOnRelPoseBetween(nullptr);
  BOOST_CHECK(!joints.first);
  BOOST_CHECK_EQUAL(joints.second->index(), ee1->index());

  function = RelativePosition::create("RelativePosition", device, ee1,
                                      JointPtr_t(), tf1, tf2);
  joints = function->dependsOnRelPoseBetween(nullptr);
  BOOST_CHECK(!joints.first);
  BOOST_CHECK_EQUAL(joints.second->index(), ee1->index());

  function = RelativeTransformationR3xSO3::create(
      "RelativeTransformationR3xSO3", device, ee1, JointPtr_t(), tf1, tf2);
  joints = function->dependsOnRelPoseBetween(nullptr);
  BOOST_CHECK(!joints.first);
  BOOST_CHECK_EQUAL(joints.second->index(), ee1->index());
  // constraint fully constrains the relative pose
  constraint = Implicit::create(
      function, ComparisonTypes_t(function->outputSpace()->nv(), EqualToZero),
      std::vector<bool>(function->outputSpace()->nv(), true));
  jointsConstrained = constraint->doesConstrainRelPoseBetween(device);
  BOOST_CHECK(!jointsConstrained.first);
  BOOST_CHECK_EQUAL(jointsConstrained.second->index(), ee1->index());

  /// test the locked joint constraint as well
  constraint = LockedJoint::create(ee2, ee2->configurationSpace()->neutral());

  std::cout << constraint->functionPtr()->name() << std::endl;
  joints = constraint->functionPtr()->dependsOnRelPoseBetween(device);
  BOOST_CHECK_EQUAL(joints.first->index(), ee2->parentJoint()->index());
  BOOST_CHECK_EQUAL(joints.second->index(), ee2->index());
  // constraint fully locks the joint
  jointsConstrained = constraint->doesConstrainRelPoseBetween(device);
  BOOST_CHECK_EQUAL(jointsConstrained.second->index(), ee2->index());
}
