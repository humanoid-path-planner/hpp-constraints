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
    config.head(offset) = se3::randomConfiguration(robot_->model());

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
  JointPtr_t ee1 = device->getJointByName ("lleg5_joint"),
             ee2 = device->getJointByName ("rleg5_joint");
  Configuration_t goal = device->currentConfiguration ();
  BOOST_REQUIRE (device);
  BasicConfigurationShooter cs (device);

  device->currentConfiguration (*cs.shoot ());
  device->computeForwardKinematics ();
  Transform3f tf1 (ee1->currentTransformation ());
  Transform3f tf2 (ee2->currentTransformation ());

  std::cout << *Orientation::create            ("Orientation"           , device, ee2, tf2)           << std::endl;
  std::cout << *Position::create               ("Position"              , device, ee2, tf2, tf1)      << std::endl;
  std::cout << *Transformation::create         ("Transformation"        , device, ee1, tf1)           << std::endl;
  std::cout << *RelativeOrientation::create    ("RelativeOrientation"   , device, ee1, ee2, tf1)      << std::endl;
  std::cout << *RelativePosition::create       ("RelativePosition"      , device, ee1, ee2, tf1, tf2) << std::endl;
  std::cout << *RelativeTransformation::create ("RelativeTransformation", device, ee1, ee2, tf1, tf2) << std::endl;
}
