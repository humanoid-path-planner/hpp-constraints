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

#include <hpp/model/configuration.hh>
#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include <hpp/model/object-factory.hh>

#include "hpp/constraints/tools.hh"

#define BOOST_TEST_MODULE StaticStability
#include <math.h>
#include <stdlib.h>

#include <boost/test/included/unit_test.hpp>
#include <limits>

using hpp::model::BodyPtr_t;
using hpp::model::Configuration_t;
using hpp::model::ConfigurationPtr_t;
using hpp::model::Device;
using hpp::model::DevicePtr_t;
using hpp::model::JointPtr_t;
using hpp::model::JointVector_t;

using std::numeric_limits;

using namespace hpp::constraints;

const static size_t NUMBER_JACOBIAN_CALCULUS = 5;
const static double HESSIAN_MAXIMUM_COEF = 1e1;
const static double DQ_MAX = 1e-2;
const static size_t MAX_NB_ERROR = 5;

static matrix3_t identity() {
  matrix3_t R;
  R.setIdentity();
  return R;
}
static fcl::Transform3f transform3f_id() {
  fcl::Transform3f T;
  T.setIdentity();
  return T;
}

hpp::model::ObjectFactory objectFactory;

class BasicConfigurationShooter {
 public:
  BasicConfigurationShooter(const DevicePtr_t& robot) : robot_(robot) {}
  virtual ConfigurationPtr_t shoot() const {
    JointVector_t jv = robot_->getJointVector();
    ConfigurationPtr_t config(new Configuration_t(robot_->configSize()));
    for (JointVector_t::const_iterator itJoint = jv.begin();
         itJoint != jv.end(); itJoint++) {
      std::size_t rank = (*itJoint)->rankInConfiguration();
      (*itJoint)->configuration()->uniformlySample(rank, *config);
    }
    return config;
  }

 private:
  const DevicePtr_t& robot_;
};  // class BasicConfigurationShooter

JointPtr_t createFreeflyerJoint(DevicePtr_t robot) {
  const std::string& name = "";
  JointPtr_t joint, parent;
  std::string jointName = name + "_xyz";
  // Translation along xyz
  joint = objectFactory.createJointTranslation3(transform3f_id());
  joint->name(jointName);
  for (size_type i = 0; i < 3; ++i) {
    joint->isBounded(i, true);
    joint->lowerBound(i, -4);
    joint->upperBound(i, +4);
  }
  robot->rootJoint(joint);
  parent = joint;

  // joint SO3
  joint = objectFactory.createJointSO3(transform3f_id());
  jointName = name + "_SO3";
  joint->name(jointName);
  parent->addChildJoint(joint);
  return joint;
}

JointPtr_t createRotationJoint(DevicePtr_t robot) {
  const std::string& name = "";
  JointPtr_t joint;
  std::string jointName = name + "_rz";

  fcl::Transform3f mat;
  mat.setIdentity();
  fcl::Matrix3f orient;
  orient(0, 0) = 0;
  orient(0, 1) = 1;
  orient(0, 2) = 0;
  orient(1, 0) = 0;
  orient(1, 1) = 0;
  orient(1, 2) = 1;
  orient(2, 0) = 1;
  orient(2, 1) = 0;
  orient(2, 2) = 0;
  mat.setRotation(orient);

  // joint rz
  JointPtr_t anchor = objectFactory.createJointAnchor(transform3f_id());
  anchor->name("anchor");
  robot->rootJoint(anchor);
  joint = objectFactory.createUnBoundedJointRotation(mat);
  joint->name(jointName);
  anchor->addChildJoint(joint);
  return joint;
}

DevicePtr_t createRobot() {
  DevicePtr_t robot = Device::create("StaticStability");
  JointPtr_t ff = createRotationJoint(robot);

  BodyPtr_t body;

  fcl::Transform3f mat;
  mat.setIdentity();
  fcl::Matrix3f orient;
  orient(0, 0) = 1;
  orient(0, 1) = 0;
  orient(0, 2) = 0;
  orient(1, 0) = 0;
  orient(1, 1) = 1;
  orient(1, 2) = 0;
  orient(2, 0) = 0;
  orient(2, 1) = 0;
  orient(2, 2) = 1;
  mat.setRotation(orient);

  // Translation along x
  JointPtr_t joint = objectFactory.createJointTranslation(mat);
  joint->name("slider");
  ff->addChildJoint(joint);
  joint->isBounded(0, true);
  joint->lowerBound(0, -4);
  joint->upperBound(0, +4);
  body = objectFactory.createBody();
  body->name("slider");
  body->mass(1);
  joint->setLinkedBody(body);
  return robot;
}

StaticStabilityPtr_t createStaticStability(DevicePtr_t d, JointPtr_t j) {
  StaticStability::Contacts_t cs;
  StaticStability::Contact_t c;

  c.joint1 = NULL;
  c.point1 = vector3_t(-1, 0, 0);
  c.normal1 = vector3_t(0, 0, 1);
  c.joint2 = j;
  c.point2 = vector3_t(-1, 0, 0);
  c.normal2 = vector3_t(0, 0, 1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(+1, 0, 0);
  c.normal1 = vector3_t(0, 0, 1);
  c.joint2 = j;
  c.point2 = vector3_t(+1, 0, 0);
  c.normal2 = vector3_t(0, 0, 1);
  cs.push_back(c);

  CenterOfMassComputationPtr_t com = CenterOfMassComputation::create(d);
  com->add(d->rootJoint());
  com->computeMass();

  StaticStabilityPtr_t fptr = StaticStability::create("test", d, cs, com);
  return fptr;
}

StaticStabilityPtr_t createStaticStabilityHard(DevicePtr_t d, JointPtr_t j) {
  StaticStability::Contacts_t cs;
  StaticStability::Contact_t c;

  c.joint1 = NULL;
  c.point1 = vector3_t(0, 0, 0);
  c.normal1 = vector3_t(0, 0, 1);
  c.joint2 = j;
  c.point2 = vector3_t(0, 0, 0);
  c.normal2 = vector3_t(0, 0, 1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(+1, 1, 0);
  c.normal1 = vector3_t(-1, 0, 0);
  c.joint2 = j;
  c.point2 = vector3_t(+1, 1, 0);
  c.normal2 = vector3_t(-1, 0, 0);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(-1, 1, 0);
  c.normal1 = vector3_t(1, 0, 0);
  c.joint2 = j;
  c.point2 = vector3_t(-1, 1, 0);
  c.normal2 = vector3_t(1, 0, 0);
  cs.push_back(c);

  CenterOfMassComputationPtr_t com = CenterOfMassComputation::create(d);
  com->add(d->rootJoint());
  com->computeMass();

  StaticStabilityPtr_t fptr = StaticStability::create("test", d, cs, com);
  return fptr;
}

StaticStabilityPtr_t createStaticStabilityHard2(DevicePtr_t d, JointPtr_t j) {
  StaticStability::Contacts_t cs;
  StaticStability::Contact_t c;

  c.joint1 = NULL;
  c.point1 = vector3_t(-1, 0, 1);
  c.normal1 = vector3_t(0, 0, -1);
  c.joint2 = j;
  c.point2 = vector3_t(-1, 0, 1);
  c.normal2 = vector3_t(0, 0, -1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(+1, 0, 1);
  c.normal1 = vector3_t(0, 0, -1);
  c.joint2 = j;
  c.point2 = vector3_t(+1, 0, 1);
  c.normal2 = vector3_t(0, 0, -1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(-1, 0, 0);
  c.normal1 = vector3_t(0, 0, 1);
  c.joint2 = j;
  c.point2 = vector3_t(-1, 0, 0);
  c.normal2 = vector3_t(0, 0, 1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(+1, 0, 0);
  c.normal1 = vector3_t(0, 0, 1);
  c.joint2 = j;
  c.point2 = vector3_t(+1, 0, 0);
  c.normal2 = vector3_t(0, 0, 1);
  cs.push_back(c);

  CenterOfMassComputationPtr_t com = CenterOfMassComputation::create(d);
  com->add(d->rootJoint());
  com->computeMass();

  StaticStabilityPtr_t fptr = StaticStability::create("test", d, cs, com);
  return fptr;
}

StaticStabilityPtr_t createStaticStabilityHard3(DevicePtr_t d, JointPtr_t j) {
  StaticStability::Contacts_t cs;
  StaticStability::Contact_t c;

  c.joint1 = NULL;
  c.point1 = vector3_t(-1, 0, 1);
  c.normal1 = vector3_t(0, 0, -1);
  c.joint2 = j;
  c.point2 = vector3_t(-1, 0, 1);
  c.normal2 = vector3_t(0, 0, -1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(0, 0, 1);
  c.normal1 = vector3_t(0, 0, -1);
  c.joint2 = j;
  c.point2 = vector3_t(0, 0, 1);
  c.normal2 = vector3_t(0, 0, -1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(+1, 0, 1);
  c.normal1 = vector3_t(0, 0, -1);
  c.joint2 = j;
  c.point2 = vector3_t(+1, 0, 1);
  c.normal2 = vector3_t(0, 0, -1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(-1, 0, 0);
  c.normal1 = vector3_t(0, 0, 1);
  c.joint2 = j;
  c.point2 = vector3_t(-1, 0, 0);
  c.normal2 = vector3_t(0, 0, 1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(+1, 0, 0);
  c.normal1 = vector3_t(0, 0, 1);
  c.joint2 = j;
  c.point2 = vector3_t(+1, 0, 0);
  c.normal2 = vector3_t(0, 0, 1);
  cs.push_back(c);

  CenterOfMassComputationPtr_t com = CenterOfMassComputation::create(d);
  com->add(d->rootJoint());
  com->computeMass();

  StaticStabilityPtr_t fptr = StaticStability::create("test", d, cs, com);
  return fptr;
}

StaticStabilityPtr_t createStaticStabilityHard4(DevicePtr_t d, JointPtr_t j) {
  StaticStability::Contacts_t cs;
  StaticStability::Contact_t c;

  c.joint1 = NULL;
  c.point1 = vector3_t(-1, 1, 1);
  c.normal1 = vector3_t(0, 0, -1);
  c.joint2 = j;
  c.point2 = vector3_t(-1, 1, 1);
  c.normal2 = vector3_t(0, 0, -1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(+1, 1, 1);
  c.normal1 = vector3_t(0, 0, -1);
  c.joint2 = j;
  c.point2 = vector3_t(+1, 1, 1);
  c.normal2 = vector3_t(0, 0, -1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(-1, 1, 0);
  c.normal1 = vector3_t(0, 0, 1);
  c.joint2 = j;
  c.point2 = vector3_t(-1, 1, 0);
  c.normal2 = vector3_t(0, 0, 1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(+1, 1, 0);
  c.normal1 = vector3_t(0, 0, 1);
  c.joint2 = j;
  c.point2 = vector3_t(+1, 1, 0);
  c.normal2 = vector3_t(0, 0, 1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(-1, 0, 1);
  c.normal1 = vector3_t(0, 0, -1);
  c.joint2 = j;
  c.point2 = vector3_t(-1, 0, 1);
  c.normal2 = vector3_t(0, 0, -1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(+1, 0, 1);
  c.normal1 = vector3_t(0, 0, -1);
  c.joint2 = j;
  c.point2 = vector3_t(+1, 0, 1);
  c.normal2 = vector3_t(0, 0, -1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(-1, 0, 0);
  c.normal1 = vector3_t(0, 0, 1);
  c.joint2 = j;
  c.point2 = vector3_t(-1, 0, 0);
  c.normal2 = vector3_t(0, 0, 1);
  cs.push_back(c);

  c.joint1 = NULL;
  c.point1 = vector3_t(+1, 0, 0);
  c.normal1 = vector3_t(0, 0, 1);
  c.joint2 = j;
  c.point2 = vector3_t(+1, 0, 0);
  c.normal2 = vector3_t(0, 0, 1);
  cs.push_back(c);

  CenterOfMassComputationPtr_t com = CenterOfMassComputation::create(d);
  com->add(d->rootJoint());
  com->computeMass();

  StaticStabilityPtr_t fptr = StaticStability::create("test", d, cs, com);
  return fptr;
}

BOOST_AUTO_TEST_CASE(static_stability) {
  DevicePtr_t device = createRobot();
  BOOST_REQUIRE(device);
  JointPtr_t rot = device->getJointByName("_rz");
  JointPtr_t slider = device->getJointByName("slider");

  Configuration_t c(3);
  c << 1, 0, 0;
  device->currentConfiguration(c);
  device->computeForwardKinematics(pinocchio::JOINT_POSITION);
  // BOOST_TEST_MESSAGE ("\"rot\" initial transform:\n" <<
  // rot->currentTransformation ()); BOOST_TEST_MESSAGE ("\"rot\" current
  // transform:\n" << rot->currentTransformation ()); BOOST_TEST_MESSAGE
  // ("\"slider\" initial transform:\n" << slider->currentTransformation ());
  // BOOST_TEST_MESSAGE ("\"slider\" current transform:\n" <<
  // slider->currentTransformation ());

  BOOST_CHECK_MESSAGE(slider->currentTransformation().isIdentity(),
                      "This transform shoud be identity:\n"
                          << slider->currentTransformation());

  StaticStabilityPtr_t fptr = createStaticStability(device, rot);
  StaticStabilityPtr_t fptrH = createStaticStabilityHard4(device, rot);
  StaticStability& f(*fptr);
  StaticStability& fH(*fptrH);
  const std::size_t nbC = 2;
  const std::size_t nbCH = 8;
  vector_t value(f.outputSize());
  vector_t valueH(fH.outputSize());
  matrix_t j(f.outputSize(), f.inputDerivativeSize());
  std::list<Configuration_t> valid, invalid;

  valid.push_back((Configuration_t(3) << 1, 0, 0).finished());
  valid.push_back((Configuration_t(3) << 1, 0, 0.5).finished());
  valid.push_back((Configuration_t(3) << 1, 0, -0.5).finished());
  valid.push_back((Configuration_t(3) << 1, 0, 1).finished());
  valid.push_back((Configuration_t(3) << 1, 0, -1).finished());
  for (std::list<Configuration_t>::const_iterator it = valid.begin();
       it != valid.end(); ++it) {
    BOOST_TEST_MESSAGE("Config " << it->transpose());
    f(value, *it);
    BOOST_TEST_MESSAGE("\"slider\" current transform:\n"
                       << slider->currentTransformation());
    BOOST_CHECK_MESSAGE(value.segment<6>(nbC).isZero(),
                        "(I - phi * phi^+) * G =\n"
                            << value.segment<6>(nbC).transpose());
    BOOST_CHECK_MESSAGE((value.segment<nbCH>(0).array() >
                         -Eigen::NumTraits<value_type>::dummy_precision())
                            .all(),
                        "Contact forces =\n"
                            << value.segment<nbC>(0).transpose());
    fH(valueH, *it);
    BOOST_CHECK_MESSAGE(valueH.segment<6>(nbCH).isZero(),
                        "(I - phi * phi^+) * G =\n"
                            << valueH.segment<6>(nbCH).transpose());
    BOOST_CHECK_MESSAGE((valueH.segment<nbCH>(0).array() >
                         -Eigen::NumTraits<value_type>::dummy_precision())
                            .all(),
                        "Contact forces =\n"
                            << valueH.segment<nbCH>(0).transpose());
  }

  BOOST_TEST_MESSAGE("Starting invalid");
  invalid.push_back((Configuration_t(3) << 0, 1, 0.5).finished());
  invalid.push_back((Configuration_t(3) << 1, 0, 2.5).finished());
  invalid.push_back((Configuration_t(3) << 1, 0, -1.5).finished());
  for (std::list<Configuration_t>::const_iterator it = invalid.begin();
       it != invalid.end(); ++it) {
    BOOST_TEST_MESSAGE("Config " << it->transpose());
    f(value, *it);
    BOOST_TEST_MESSAGE("\"slider\" current transform:\n"
                       << slider->currentTransformation());
    BOOST_CHECK_MESSAGE(
        !value.segment<6>(nbC).isZero() ||
            (value.segment<nbC>(0).array() <
             -Eigen::NumTraits<value_type>::dummy_precision())
                .any(),
        "Should not be stable:\n"
        "(I - phi * phi^+) * G =\n"
            << value.segment<6>(nbC).transpose() << "\nContact forces =\n"
            << value.segment<nbC>(0).transpose() << "\nJoint transform:\n"
            << slider->currentTransformation());

    fH(valueH, *it);
    BOOST_CHECK_MESSAGE(
        !valueH.segment<6>(nbCH).isZero() ||
            (valueH.segment<nbCH>(0).array() <
             -Eigen::NumTraits<value_type>::dummy_precision())
                .any(),
        "Should not be stable:\n"
        "(I - phi * phi^+) * G =\n"
            << valueH.segment<6>(nbCH).transpose() << "\nContact forces =\n"
            << valueH.segment<nbCH>(0).transpose() << "\nJoint transform:\n"
            << slider->currentTransformation());
  }
}

BOOST_AUTO_TEST_CASE(static_stability_phi) {
  DevicePtr_t device = createRobot();
  BOOST_REQUIRE(device);
  JointPtr_t rot = device->getJointByName("_rz");
  JointPtr_t slider = device->getJointByName("slider");

  Configuration_t c(3);
  c << 1, 0, 0;
  device->currentConfiguration(c);
  device->computeForwardKinematics(pinocchio::JOINT_POSITION);

  BOOST_CHECK_MESSAGE(slider->currentTransformation().isIdentity(),
                      "This transform shoud be identity:\n"
                          << slider->currentTransformation());

  StaticStabilityPtr_t fptr =
      createStaticStabilityHard4(device, device->rootJoint());
  StaticStability& f(*fptr);
  const std::size_t nbC = 8;
  vector_t value(f.outputSize());
  matrix_t j(f.outputSize(), f.inputDerivativeSize());
  std::list<Configuration_t> valid, invalid;

  const double sqr2 = sqrt(2);
  const double invsr2 = 1 / sqr2;
  Configuration_t center(3);
  center << 1, 0, 0;
  vector_t FCenter(nbC);
  FCenter << 0, 0, 0, 0, 0, 0, 1, 1;
  FCenter /= 2;
  Configuration_t point2(3);
  point2 << -invsr2, invsr2, sqr2;
  vector_t F2(nbC);
  F2.setZero();
  F2[2] = 1;
  Configuration_t point3(3);
  point3 << invsr2, invsr2, sqr2;
  vector_t F3(nbC);
  F3.setZero();
  F3[3] = 1;
  Configuration_t point6(3);
  point6 << -1, 0, 1;
  vector_t F6(nbC);
  F6.setZero();
  F6[6] = 1;
  Configuration_t point7(3);
  point7 << 1, 0, 1;
  vector_t F7(nbC);
  F7.setZero();
  F7[7] = 1;

  f(value, center);
  BOOST_CHECK_MESSAGE((value.segment(0, nbC).array() >= 0).all(),
                      "No positive solution found:\n"
                          << value);
  BOOST_CHECK_MESSAGE(value.segment<6>(nbC).isZero(), "No solution found:\n"
                                                          << value);
  BOOST_CHECK_MESSAGE(
      (f.phi().value() * FCenter + StaticStability::Gravity).isZero(),
      "Residual is:\n"
          << f.phi().value() * FCenter);
  f(value, point2);
  BOOST_CHECK_MESSAGE((value.segment(0, nbC).array() >= 0).all(),
                      "No positive solution found:\n"
                          << value);
  BOOST_CHECK_MESSAGE(value.segment<6>(nbC).isZero(), "No solution found:\n"
                                                          << value);
  BOOST_CHECK_MESSAGE(
      (f.phi().value() * F2 + StaticStability::Gravity).isZero(),
      "Residual is:\n"
          << f.phi().value() * F2);
  f(value, point3);
  BOOST_CHECK_MESSAGE((value.segment(0, nbC).array() >= 0).all(),
                      "No positive solution found:\n"
                          << value);
  BOOST_CHECK_MESSAGE(value.segment<6>(nbC).isZero(), "No solution found:\n"
                                                          << value);
  BOOST_CHECK_MESSAGE(
      (f.phi().value() * F3 + StaticStability::Gravity).isZero(),
      "Residual is:\n"
          << f.phi().value() * F3);
  f(value, point6);
  BOOST_CHECK_MESSAGE((value.segment(0, nbC).array() >= 0).all(),
                      "No positive solution found:\n"
                          << value);
  BOOST_CHECK_MESSAGE(value.segment<6>(nbC).isZero(), "No solution found:\n"
                                                          << value);
  BOOST_CHECK_MESSAGE(
      (f.phi().value() * F6 + StaticStability::Gravity).isZero(),
      "Residual is:\n"
          << f.phi().value() * F6);
  f(value, point7);
  BOOST_CHECK_MESSAGE((value.segment(0, nbC).array() >= 0).all(),
                      "No positive solution found:\n"
                          << value);
  BOOST_CHECK_MESSAGE(value.segment<6>(nbC).isZero(), "No solution found:\n"
                                                          << value);
  BOOST_CHECK_MESSAGE(
      (f.phi().value() * F7 + StaticStability::Gravity).isZero(),
      "Residual is:\n"
          << f.phi().value() * F7);
}
