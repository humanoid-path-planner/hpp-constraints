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

#define EIGEN_RUNTIME_NO_MALLOC

#include "hpp/constraints/generic-transformation.hh"

#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/simple-device.hh>

#include "hpp/constraints/tools.hh"

#define BOOST_TEST_MODULE hpp_constraints
#include <boost/test/included/unit_test.hpp>

#include <stdlib.h>

using hpp::pinocchio::Configuration_t;
using hpp::pinocchio::ConfigurationPtr_t;
using hpp::pinocchio::Device;
using hpp::pinocchio::DevicePtr_t;
using hpp::pinocchio::JointPtr_t;
using hpp::pinocchio::Transform3f;

using namespace hpp::constraints;

class BasicConfigurationShooter
{
public:
  BasicConfigurationShooter (const DevicePtr_t& robot) : robot_ (robot)
  {
  }
  virtual ConfigurationPtr_t shoot () const
  {
    size_type extraDim = robot_->extraConfigSpace ().dimension ();
    size_type offset = robot_->configSize () - extraDim;

    Configuration_t config(robot_->configSize ());
    config.head(offset) = ::pinocchio::randomConfiguration(robot_->model());

    // Shoot extra configuration variables
    for (size_type i=0; i<extraDim; ++i) {
      value_type lower = robot_->extraConfigSpace ().lower (i);
      value_type upper = robot_->extraConfigSpace ().upper (i);
      value_type range = upper - lower;
      if ((range < 0) ||
          (range == std::numeric_limits<double>::infinity())) {
        std::ostringstream oss
          ("Cannot uniformy sample extra config variable ");
        oss << i << ". min = " <<lower<< ", max = " << upper << std::endl;
        throw std::runtime_error (oss.str ());
      }
      config [offset + i] = lower + (upper - lower) * rand ()/RAND_MAX;
    }
    return boost::make_shared<Configuration_t>(config);
  }
private:
  const DevicePtr_t& robot_;
}; // class BasicConfigurationShooter

BOOST_AUTO_TEST_CASE (print) {
  DevicePtr_t device = hpp::pinocchio::humanoidSimple ("test");
  device->model().upperPositionLimit.head<3>().setOnes();
  device->model().lowerPositionLimit.head<3>().setZero();
  JointPtr_t ee1 = device->getJointByName ("lleg5_joint"),
             ee2 = device->getJointByName ("rleg5_joint");
  BOOST_REQUIRE (device);
  BasicConfigurationShooter cs (device);

  device->currentConfiguration (*cs.shoot ());
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

  Configuration_t q1 = *cs.shoot(),
                  q2 = *cs.shoot();
  for (std::size_t i = 0; i < functions.size(); ++i) {
    DifferentiableFunctionPtr_t f = functions[i];

    std::cout << *f << std::endl;

    LiegroupElement v (f->outputSpace());
    matrix_t J (f->outputDerivativeSize(), f->inputDerivativeSize());

    f->value    (v, q1);
    f->jacobian (J, q1);
    // TODO this is broken at the moment because of the introduction
    // of a multithreaded device.
    // Eigen::internal::set_is_malloc_allowed(false);
    f->value    (v, q2);
    f->jacobian (J, q2);
    // Eigen::internal::set_is_malloc_allowed(true);
  }

  // Check active parameters
  ArrayXb ap1 = Orientation::create ("Orientation"           , device, ee1, tf1)->activeParameters();
  ArrayXb ap2 = Orientation::create ("Orientation"           , device, ee2, tf2)->activeParameters();
  ArrayXb ap12 = RelativeOrientation::create    ("RelativeOrientation"   , device, ee1, ee2, tf1)->activeParameters();

  ArrayXb not_ap1 = (ap1 == false);
  ArrayXb not_ap2 = (ap2 == false);
  BOOST_CHECK ((ap12 == ((not_ap1 && ap2) || (ap1 && not_ap2))).all());

  // Check active derivative parameters
  ap1 = Orientation::create ("Orientation"           , device, ee1, tf1)->activeDerivativeParameters();
  ap2 = Orientation::create ("Orientation"           , device, ee2, tf2)->activeDerivativeParameters();
  ap12 = RelativeOrientation::create    ("RelativeOrientation"   , device, ee1, ee2, tf1)->activeDerivativeParameters();

  not_ap1 = (ap1 == false);
  not_ap2 = (ap2 == false);
  BOOST_CHECK ((ap12 == ((not_ap1 && ap2) || (ap1 && not_ap2))).all());
}

BOOST_AUTO_TEST_CASE (multithread) {
  DevicePtr_t device = hpp::pinocchio::humanoidSimple ("test");
  device->numberDeviceData (4);
  device->model().upperPositionLimit.head<3>().setOnes();
  device->model().lowerPositionLimit.head<3>().setZero();
  JointPtr_t ee1 = device->getJointByName ("lleg5_joint"),
             ee2 = device->getJointByName ("rleg5_joint");
  BOOST_REQUIRE (device);
  BasicConfigurationShooter cs (device);

  device->currentConfiguration (*cs.shoot ());
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

  const int N = 10;
  Configuration_t q = *cs.shoot();
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
