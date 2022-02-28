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

#define BOOST_TEST_MODULE multithread_test

#include <pinocchio/fwd.hpp>

#include <boost/test/included/unit_test.hpp>

#include <stdlib.h>

#include <hpp/constraints/configuration-constraint.hh>
#include <hpp/constraints/convex-shape-contact.hh>
#include <hpp/constraints/generic-transformation.hh>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/simple-device.hh>

#include <../tests/util.hh>
#include <../tests/convex-shape-contact-function.hh>

using hpp::pinocchio::Configuration_t;
using hpp::pinocchio::ConfigurationPtr_t;
using hpp::pinocchio::Device;
using hpp::pinocchio::DevicePtr_t;
using hpp::pinocchio::JointPtr_t;
using hpp::pinocchio::Transform3f;

using namespace hpp::constraints;

BOOST_AUTO_TEST_CASE (multithread) {
  DevicePtr_t device = hpp::pinocchio::unittest::makeDevice(
      hpp::pinocchio::unittest::HumanoidSimple);
  device->numberDeviceData (4);
  JointPtr_t ee1 = device->getJointByName ("lleg5_joint"),
             ee2 = device->getJointByName ("rleg5_joint");
  BOOST_REQUIRE (device);

  Configuration_t q;
  randomConfig (device, q);
  device->currentConfiguration (q);
  device->computeForwardKinematics ();
  Transform3f tf1 (ee1->currentTransformation ());
  Transform3f tf2 (ee2->currentTransformation ());

  std::vector<DifferentiableFunctionPtr_t> functions;
  functions.push_back(ConfigurationConstraint::create ("Configuration", device, device->currentConfiguration()));
  functions.push_back(Orientation::create            ("Orientation"           , device, ee2, tf2)          );
  functions.push_back(Position::create               ("Position"              , device, ee2, tf2, tf1)     );
  functions.push_back(Transformation::create         ("Transformation"        , device, ee1, tf1)          );
  functions.push_back(RelativeOrientation::create    ("RelativeOrientation"   , device, ee1, ee2, tf1)     );
  functions.push_back(RelativePosition::create       ("RelativePosition"      , device, ee1, ee2, tf1, tf2));
  functions.push_back(RelativeTransformation::create ("RelativeTransformation", device, ee1, ee2, tf1, tf2));
  functions.push_back(createConvexShapeContact_triangles (device, ee1, "ConvexShapeContact triangle"));
  functions.push_back(createConvexShapeContact_punctual  (device, ee1, "ConvexShapeContact punctual"));
  functions.push_back(createConvexShapeContact_convex    (device, ee1, "ConvexShapeContact convex"));

  const int N = 100;
  randomConfig (device, q);
  for (std::size_t i = 0; i < functions.size(); ++i) {
    DifferentiableFunctionPtr_t f = functions[i];

    std::vector <LiegroupElement> vs (N, LiegroupElement (f->outputSpace()));
    std::vector <matrix_t> Js (N, matrix_t(f->outputDerivativeSize(), f->inputDerivativeSize()));
#pragma omp parallel for
    for (int j = 0; j < N; ++j) {
      f->value    (vs[j], q);
      f->jacobian (Js[j], q);
    }

    for (int j = 1; j < N; ++j) {
      BOOST_CHECK_EQUAL (vs[0].vector(), vs[j].vector());
      BOOST_CHECK_EQUAL (Js[0]         , Js[j]);
    }
  }
}
