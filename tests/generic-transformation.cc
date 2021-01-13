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
#include <hpp/constraints/solver/by-substitution.hh>

#include <sstream>
#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/joint-collection.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/simple-device.hh>
#include <hpp/pinocchio/serialization.hh>
#include <hpp/pinocchio/urdf/util.hh>

#include "hpp/constraints/tools.hh"

#define BOOST_TEST_MODULE hpp_constraints
#include <boost/test/included/unit_test.hpp>

#include <stdlib.h>

using hpp::pinocchio::Configuration_t;
using hpp::pinocchio::ConfigurationPtr_t;
using hpp::pinocchio::Device;
using hpp::pinocchio::DevicePtr_t;
using hpp::pinocchio::JointPtr_t;
using hpp::pinocchio::LiegroupSpace;
using hpp::pinocchio::Transform3f;

using hpp::pinocchio::urdf::loadModelFromString;

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
  DevicePtr_t robot_;
}; // class BasicConfigurationShooter

BOOST_AUTO_TEST_CASE (print) {
  DevicePtr_t device = hpp::pinocchio::unittest::makeDevice(
      hpp::pinocchio::unittest::HumanoidSimple);
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
  DevicePtr_t device = hpp::pinocchio::unittest::makeDevice(
      hpp::pinocchio::unittest::HumanoidSimple);
  device->numberDeviceData (4);
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
  functions.push_back(RelativeOrientation::create    ("RelativeOrientation"   , device, ee1, JointPtr_t(), tf1)     );
  functions.push_back(RelativePosition::create       ("RelativePosition"      , device, ee1, JointPtr_t(), tf1, tf2));
  functions.push_back(RelativeTransformation::create ("RelativeTransformation", device, ee1, JointPtr_t(), tf1, tf2));

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

BOOST_AUTO_TEST_CASE (serialization) {
  DevicePtr_t device = hpp::pinocchio::unittest::makeDevice(
      hpp::pinocchio::unittest::HumanoidSimple);
  device->numberDeviceData (4);
  JointPtr_t ee1 = device->getJointByName ("lleg5_joint"),
             ee2 = device->getJointByName ("rleg5_joint");
  BOOST_REQUIRE (device);

  device->currentConfiguration (device->neutralConfiguration());
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
  functions.push_back(RelativeOrientation::create    ("RelativeOrientation"   , device, ee1, JointPtr_t(), tf1)     );
  functions.push_back(RelativePosition::create       ("RelativePosition"      , device, ee1, JointPtr_t(), tf1, tf2));
  functions.push_back(RelativeTransformation::create ("RelativeTransformation", device, ee1, JointPtr_t(), tf1, tf2));

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

BOOST_AUTO_TEST_CASE(RelativeTransformation_SE3)
{
  const std::string model("<robot name=\"box\">"
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
  for (std::size_t i=0; i<2; ++i) {
    vector_t l(7); l << -2,-2,-2,-1,-1,-1,-1;
    vector_t u(7); u <<  2, 2, 2, 1, 1, 1, 1;
    robot->jointAt(i)->lowerBounds(l);
    robot->jointAt(i)->upperBounds(u);
  }
  // Create constraint
  //
  // Frame in joint 1
  //   R =    0.707107, -0.707107,         0
  //          0,         0,                1
  //         -0.707107, -0.707107,         0
  //   p =         0.1,       0,       -0.03
  //
  // Frame in joint 2
  //   R =   1, 0, 0
  //         0, 1, 0
  //         0, 0, 1
  //   p =    0,    0,    -0.2

  matrix3_t R1, R2;
  vector3_t p1, p2;
  value_type a(sqrt(2)/2);
  R1 << a, -a, 0,
    0, 0, 1,
    -a, -a, 0;
  p1 << 0.1, 0, -0.03;
  R2.setIdentity();
  p2 << 0, 0, -0.2;
  Transform3f tf1(R1, p1), tf2(R2,p2);
  std::vector<bool> mask = {false, false, false, false, false, true};
  ImplicitPtr_t constraint
    (Implicit::create (RelativeTransformationSE3::create
		       ("RelativeTransformationSE3", robot, j1, j2,
			tf1, tf2),
		       6 * Equality, mask));
  BasicConfigurationShooter cs (robot);
  solver::BySubstitution solver(robot->configSpace());
  solver.errorThreshold(1e-10);
  solver.add(constraint);
  for (size_type i=0; i<10; ++i)
  {
    ConfigurationPtr_t q(cs.shoot());
    LiegroupElement h(LiegroupSpace::R3xSO3());
    vector6_t error;
    solver.rightHandSideFromConfig(*q);
    BOOST_CHECK(solver.isSatisfied(*q, error));
    constraint->function().value(h, *q);
    std::cout << "q =  " << q->transpose() << std::endl;
    std::cout << "h(q) = " << h << std::endl;
    std::cout << "error = " << error.transpose() << std::endl;
    std::cout << solver << std::endl;
  }
}
