// Copyright (c) 2023, CNRS
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

#include <hpp/constraints/generic-transformation.hh>
#include <hpp/constraints/relative-com.hh>
#include <hpp/pinocchio/simple-device.hh>

#define BOOST_TEST_MODULE RelativeCom
#include <boost/test/included/unit_test.hpp>

using hpp::pinocchio::Configuration_t;
using hpp::pinocchio::DevicePtr_t;
using hpp::pinocchio::JointPtr_t;
using hpp::pinocchio::LiegroupElement;
using hpp::pinocchio::LiegroupSpacePtr_t;
using hpp::pinocchio::Transform3s;
using hpp::pinocchio::vector3_t;
using hpp::pinocchio::unittest::HumanoidRomeo;
using hpp::pinocchio::unittest::makeDevice;

using hpp::constraints::DifferentiableFunctionPtr_t;
using hpp::constraints::GenericTransformation;
using hpp::constraints::OrientationBit;
using hpp::constraints::PositionBit;
using hpp::constraints::RelativeBit;
using hpp::constraints::RelativeCom;
using hpp::constraints::RelativeComPtr_t;

BOOST_AUTO_TEST_CASE(relative_com) {
  DevicePtr_t robot(makeDevice(HumanoidRomeo));
  JointPtr_t leftAnkle(robot->getJointByName("LAnkleRoll"));
  JointPtr_t rightAnkle(robot->getJointByName("RAnkleRoll"));
  assert(leftAnkle);
  assert(rightAnkle);
  // Create relative com constraint
  vector3_t reference;
  reference.setZero();
  RelativeComPtr_t rc(
      RelativeCom::create("relative-com", robot, leftAnkle, reference));
  Configuration_t q0(robot->neutralConfiguration());
  // Create relative pose constraint between ankles.
  Transform3s ref1, ref2;
  ref1.setIdentity();
  ref2.setIdentity();
  ref2.translation() << 0, -0.2, 0;
  std::vector<bool> mask(6, true);
  DifferentiableFunctionPtr_t fp(
      GenericTransformation<PositionBit | OrientationBit | RelativeBit>::create(
          "foot-pose", robot, leftAnkle, rightAnkle, ref1, ref2, mask));

  // Get function values in neutral configuration.
  LiegroupElement rc0(rc->outputSpace());
  LiegroupElement fp0(fp->outputSpace());

  rc->value(rc0, q0);
  fp->value(fp0, q0);

  // Translate root joint
  Configuration_t q1(q0);
  q1.head<3>() << 1., 2., 3.;
  LiegroupElement rc1(rc->outputSpace());
  LiegroupElement fp1(fp->outputSpace());
  rc->value(rc1, q1);
  fp->value(fp1, q1);
  BOOST_CHECK((rc1 - rc0).norm() < 1e-8);
  BOOST_CHECK((fp1 - fp0).norm() < 1e-8);
}
