// Copyright (c) 2016, Joseph Mirabel
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

#define BOOST_TEST_MODULE multithread_test
#include <boost/test/included/unit_test.hpp>

#include <stdlib.h>

#include <hpp/constraints/generic-transformation.hh>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
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
  DevicePtr_t device = hpp::pinocchio::humanoidSimple ("test");
  device->numberDeviceData (4);
  device->model().upperPositionLimit.head<3>().setOnes();
  device->model().lowerPositionLimit.head<3>().setZero();
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
  functions.push_back(Orientation::create            ("Orientation"           , device, ee2, tf2)          );
  functions.push_back(Position::create               ("Position"              , device, ee2, tf2, tf1)     );
  functions.push_back(Transformation::create         ("Transformation"        , device, ee1, tf1)          );
  functions.push_back(RelativeOrientation::create    ("RelativeOrientation"   , device, ee1, ee2, tf1)     );
  functions.push_back(RelativePosition::create       ("RelativePosition"      , device, ee1, ee2, tf1, tf2));
  functions.push_back(RelativeTransformation::create ("RelativeTransformation", device, ee1, ee2, tf1, tf2));
  functions.push_back(createConvexShapeContact_triangles (device, ee1, "ConvexShapeContact triangle"));
  functions.push_back(createConvexShapeContact_punctual  (device, ee1, "ConvexShapeContact punctual"));
  functions.push_back(createConvexShapeContact_convex    (device, ee1, "ConvexShapeContact convex"));

  const int N = 10;
  randomConfig (device, q);
  for (std::size_t i = 0; i < functions.size(); ++i) {
    DifferentiableFunctionPtr_t f = functions[i];

    std::vector <LiegroupElement> vs (N, LiegroupElement (f->outputSpace()));
    std::vector <matrix_t> Js (N, matrix_t(f->outputDerivativeSize(), f->inputDerivativeSize()));
#pragma omp parallel for
    for (int j = 0; j < 10; ++j) {
      f->value    (vs[j], q);
      f->jacobian (Js[j], q);
    }

    for (int j = 1; j < N; ++j) {
      BOOST_CHECK_EQUAL (vs[0].vector(), vs[j].vector());
      BOOST_CHECK_EQUAL (Js[0]         , Js[j]);
    }
  }
}
