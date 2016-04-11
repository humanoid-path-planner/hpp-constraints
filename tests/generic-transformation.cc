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

#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include <hpp/model/configuration.hh>
#include <hpp/model/object-factory.hh>

#include "hpp/constraints/position.hh"
#include "hpp/constraints/orientation.hh"
#include "hpp/constraints/transformation.hh"
#include "hpp/constraints/relative-position.hh"
#include "hpp/constraints/relative-orientation.hh"
#include "hpp/constraints/relative-transformation.hh"
#include "hpp/constraints/generic-transformation.hh"
#include "hpp/constraints/tools.hh"

#define BOOST_TEST_MODULE hpp_constraints
#include <boost/test/included/unit_test.hpp>

#include <stdlib.h>
#include <limits>
#include <math.h>

using hpp::model::Configuration_t;
using hpp::model::ConfigurationPtr_t;
using hpp::model::Device;
using hpp::model::DevicePtr_t;
using hpp::model::JointPtr_t;
using hpp::model::JointVector_t;
using hpp::model::BodyPtr_t;

using std::numeric_limits;
using boost::assign::list_of;

using namespace hpp::constraints;

const static size_t NUMBER_JACOBIAN_CALCULUS = 5;
const static double HESSIAN_MAXIMUM_COEF = 1e1;
const static double DQ_MAX = 1e-2;
const static size_t MAX_NB_ERROR = 5;

static matrix3_t identity () { matrix3_t R; R.setIdentity (); return R;}

hpp::model::ObjectFactory objectFactory;

class BasicConfigurationShooter
{
public:
  BasicConfigurationShooter (const DevicePtr_t& robot) : robot_ (robot)
  {
  }
  virtual ConfigurationPtr_t shoot () const
  {
    JointVector_t jv = robot_->getJointVector ();
    ConfigurationPtr_t config (new Configuration_t (robot_->configSize ()));
    for (JointVector_t::const_iterator itJoint = jv.begin ();
	 itJoint != jv.end (); itJoint++) {
      std::size_t rank = (*itJoint)->rankInConfiguration ();
      (*itJoint)->configuration ()->uniformlySample (rank, *config);
    }
    return config;
  }
private:
  const DevicePtr_t& robot_;
}; // class BasicConfigurationShooter

JointPtr_t createFreeflyerJoint (DevicePtr_t robot)
{
  const std::string& name = robot->name ();
  fcl::Transform3f mat; mat.setIdentity ();
  JointPtr_t joint, parent;
  const fcl::Vec3f T = mat.getTranslation ();
  std::string jointName = name + "_x";
  // Translation along x
  fcl::Matrix3f permutation;
  joint = objectFactory.createJointTranslation (mat);
  joint->name (jointName);
  joint->isBounded (0, 1);
  joint->lowerBound (0, -4);
  joint->upperBound (0, +4);
  robot->rootJoint (joint);
  parent = joint;

  // Translation along y
  permutation (0,0) = 0; permutation (0,1) = -1; permutation (0,2) = 0;
  permutation (1,0) = 1; permutation (1,1) =  0; permutation (1,2) = 0;
  permutation (2,0) = 0; permutation (2,1) =  0; permutation (2,2) = 1;
  fcl::Transform3f pos;
  pos.setRotation (permutation * mat.getRotation ());
  pos.setTranslation (T);
  joint = objectFactory.createJointTranslation (pos);
  jointName = name + "_y";
  joint->name (jointName);
  joint->isBounded (0, 1);
  joint->lowerBound (0, -4);
  joint->upperBound (0, +4);
  parent->addChildJoint (joint);
  parent = joint;

  // Translation along z
  permutation (0,0) = 0; permutation (0,1) = 0; permutation (0,2) = -1;
  permutation (1,0) = 0; permutation (1,1) = 1; permutation (1,2) =  0;
  permutation (2,0) = 1; permutation (2,1) = 0; permutation (2,2) =  0;
  pos.setRotation (permutation * mat.getRotation ());
  pos.setTranslation (T);
  joint = objectFactory.createJointTranslation (pos);
  jointName = name + "_z";
  joint->name (jointName);
  joint->isBounded (0, 1);
  joint->lowerBound (0, -4);
  joint->upperBound (0, +4);
  parent->addChildJoint (joint);
  parent = joint;
  // joint SO3
  joint = objectFactory.createJointSO3 (mat);
  jointName = name + "_SO3";
  joint->name (jointName);
  parent->addChildJoint (joint);
  return joint;
}

DevicePtr_t createRobot ()
{
  DevicePtr_t robot = Device::create ("test");
  JointPtr_t waist = createFreeflyerJoint (robot);
  JointPtr_t parent = waist;
  BodyPtr_t body;
  fcl::Transform3f pos;
  fcl::Matrix3f orient;
  // Right leg joint 0
  orient (0,0) = 0; orient (0,1) = 0; orient (0,2) = -1;
  orient (1,0) = 0; orient (1,1) = 1; orient (1,2) =  0;
  orient (2,0) = 1; orient (2,1) = 0; orient (2,2) =  0;
  pos.setRotation (orient);
  pos.setTranslation (fcl::Vec3f (0, -0.08, 0));
  JointPtr_t joint = objectFactory.createBoundedJointRotation (pos);
  joint->name ("RLEG_0");
  parent->addChildJoint (joint);
  body = objectFactory.createBody ();
  body->name ("RLEG_BODY0");
  joint->setLinkedBody (body);
  parent = joint;
  // Right leg joint 1
  orient (0,0) = 1; orient (0,1) = 0; orient (0,2) = 0;
  orient (1,0) = 0; orient (1,1) = 1; orient (1,2) = 0;
  orient (2,0) = 0; orient (2,1) = 0; orient (2,2) = 1;
  pos.setRotation (orient);
  pos.setTranslation (fcl::Vec3f (0, -0.08, 0));
  joint = objectFactory.createBoundedJointRotation (pos);
  joint->name ("RLEG_1");
  parent->addChildJoint (joint);
  body = objectFactory.createBody ();
  body->name ("RLEG_BODY1");
  joint->setLinkedBody (body);
  parent = joint;
  // Right leg joint 2
  orient (0,0) = 0; orient (0,1) = -1; orient (0,2) = 0;
  orient (1,0) = 1; orient (1,1) =  0; orient (1,2) = 0;
  orient (2,0) = 0; orient (2,1) =  0; orient (2,2) = 1;
  pos.setRotation (orient);
  pos.setTranslation (fcl::Vec3f (0, -0.08, 0));
  joint = objectFactory.createBoundedJointRotation (pos);
  joint->name ("RLEG_2");
  parent->addChildJoint (joint);
  body = objectFactory.createBody ();
  body->name ("RLEG_BODY2");
  joint->setLinkedBody (body);
  parent = joint;
  // Right leg joint 3: knee
  orient (0,0) = 0; orient (0,1) = -1; orient (0,2) = 0;
  orient (1,0) = 1; orient (1,1) =  0; orient (1,2) = 0;
  orient (2,0) = 0; orient (2,1) =  0; orient (2,2) = 1;
  pos.setRotation (orient);
  pos.setTranslation (fcl::Vec3f (0, -0.08, -0.35));
  joint = objectFactory.createBoundedJointRotation (pos);
  joint->name ("RLEG_3");
  parent->addChildJoint (joint);
  body = objectFactory.createBody ();
  body->name ("RLEG_BODY3");
  joint->setLinkedBody (body);
  parent = joint;
  // Right leg joint 4: ankle
  orient (0,0) = 0; orient (0,1) = -1; orient (0,2) = 0;
  orient (1,0) = 1; orient (1,1) =  0; orient (1,2) = 0;
  orient (2,0) = 0; orient (2,1) =  0; orient (2,2) = 1;
  pos.setRotation (orient);
  pos.setTranslation (fcl::Vec3f (0, -0.08, -0.70));
  joint = objectFactory.createBoundedJointRotation (pos);
  joint->name ("RLEG_4");
  parent->addChildJoint (joint);
  body = objectFactory.createBody ();
  body->name ("RLEG_BODY4");
  joint->setLinkedBody (body);
  parent = joint;
  // Right leg joint 5: ankle
  orient (0,0) = 1; orient (0,1) = 0; orient (0,2) = 0;
  orient (1,0) = 0; orient (1,1) = 1; orient (1,2) = 0;
  orient (2,0) = 0; orient (2,1) = 0; orient (2,2) = 1;
  pos.setRotation (orient);
  pos.setTranslation (fcl::Vec3f (0, -0.08, -0.70));
  joint = objectFactory.createBoundedJointRotation (pos);
  joint->name ("RLEG_5");
  parent->addChildJoint (joint);
  body = objectFactory.createBody ();
  body->name ("RLEG_BODY5");
  joint->setLinkedBody (body);

  // Left leg
  parent = waist;
  // Left leg joint 0
  orient (0,0) = 0; orient (0,1) = 0; orient (0,2) = -1;
  orient (1,0) = 0; orient (1,1) = 1; orient (1,2) =  0;
  orient (2,0) = 1; orient (2,1) = 0; orient (2,2) =  0;
  pos.setRotation (orient);
  pos.setTranslation (fcl::Vec3f (0, -0.08, 0));
  joint = objectFactory.createBoundedJointRotation (pos);
  joint->name ("LLEG_0");
  parent->addChildJoint (joint);
  body = objectFactory.createBody ();
  body->name ("LLEG_BODY0");
  joint->setLinkedBody (body);
  parent = joint;
  // Left leg joint 1
  orient (0,0) = 1; orient (0,1) = 0; orient (0,2) = 0;
  orient (1,0) = 0; orient (1,1) = 1; orient (1,2) = 0;
  orient (2,0) = 0; orient (2,1) = 0; orient (2,2) = 1;
  pos.setRotation (orient);
  pos.setTranslation (fcl::Vec3f (0, -0.08, 0));
  joint = objectFactory.createBoundedJointRotation (pos);
  joint->name ("LLEG_1");
  parent->addChildJoint (joint);
  body = objectFactory.createBody ();
  body->name ("LLEG_BODY1");
  joint->setLinkedBody (body);
  parent = joint;
  // Left leg joint 2
  orient (0,0) = 0; orient (0,1) = -1; orient (0,2) = 0;
  orient (1,0) = 1; orient (1,1) =  0; orient (1,2) = 0;
  orient (2,0) = 0; orient (2,1) =  0; orient (2,2) = 1;
  pos.setRotation (orient);
  pos.setTranslation (fcl::Vec3f (0, -0.08, 0));
  joint = objectFactory.createBoundedJointRotation (pos);
  joint->name ("LLEG_2");
  parent->addChildJoint (joint);
  body = objectFactory.createBody ();
  body->name ("LLEG_BODY2");
  joint->setLinkedBody (body);
  parent = joint;
  // Left leg joint 3: knee
  orient (0,0) = 0; orient (0,1) = -1; orient (0,2) = 0;
  orient (1,0) = 1; orient (1,1) =  0; orient (1,2) = 0;
  orient (2,0) = 0; orient (2,1) =  0; orient (2,2) = 1;
  pos.setRotation (orient);
  pos.setTranslation (fcl::Vec3f (0, -0.08, -0.35));
  joint = objectFactory.createBoundedJointRotation (pos);
  joint->name ("LLEG_3");
  parent->addChildJoint (joint);
  body = objectFactory.createBody ();
  body->name ("LLEG_BODY3");
  joint->setLinkedBody (body);
  parent = joint;
  // Left leg joint 4: ankle
  orient (0,0) = 0; orient (0,1) = -1; orient (0,2) = 0;
  orient (1,0) = 1; orient (1,1) =  0; orient (1,2) = 0;
  orient (2,0) = 0; orient (2,1) =  0; orient (2,2) = 1;
  pos.setRotation (orient);
  pos.setTranslation (fcl::Vec3f (0, -0.08, -0.70));
  joint = objectFactory.createBoundedJointRotation (pos);
  joint->name ("LLEG_4");
  parent->addChildJoint (joint);
  body = objectFactory.createBody ();
  body->name ("LLEG_BODY4");
  joint->setLinkedBody (body);
  parent = joint;
  // Left leg joint 5: ankle
  orient (0,0) = 1; orient (0,1) = 0; orient (0,2) = 0;
  orient (1,0) = 0; orient (1,1) = 1; orient (1,2) = 0;
  orient (2,0) = 0; orient (2,1) = 0; orient (2,2) = 1;
  pos.setRotation (orient);
  pos.setTranslation (fcl::Vec3f (0, -0.08, -0.70));
  joint = objectFactory.createBoundedJointRotation (pos);
  joint->name ("LLEG_5");
  parent->addChildJoint (joint);
  body = objectFactory.createBody ();
  body->name ("LLEG_BODY5");
  joint->setLinkedBody (body);

  return robot;
}

void timings (DevicePtr_t dev, DifferentiableFunctionPtr_t f,
    const std::size_t iter);

void check_consistent (DevicePtr_t dev,
    DifferentiableFunctionPtr_t f, DifferentiableFunctionPtr_t g,
    const value_type alpha = 1)
{
  BasicConfigurationShooter cs (dev);
  // std::cout << f->name() << '\n' << g->name() << '\n';
  BOOST_CHECK(f->outputSize()==g->outputSize());

  vector_t value1 = vector_t (f->outputSize ());
  vector_t value2 = vector_t (g->outputSize ());
  matrix_t jacobian1 = matrix_t (f->outputSize (), dev->numberDof ());
  matrix_t jacobian2 = matrix_t (g->outputSize (), dev->numberDof ());
  for (size_t i = 0; i < NUMBER_JACOBIAN_CALCULUS; i++) {
    ConfigurationPtr_t q = cs.shoot ();
    (*f) (value1, *q);
    (*g) (value2, *q);
    vector_t d = value2 - alpha * value1;
    // std::cout << d.transpose() << std::endl;
    BOOST_CHECK(value1.isApprox(alpha*value2));
    f->jacobian (jacobian1, *q);
    g->jacobian (jacobian2, *q);
    matrix_t diffJ = jacobian2 - alpha*jacobian1;
    // std::cout << diffJ.norm() << std::endl;
    BOOST_CHECK(jacobian1.isApprox(alpha*jacobian2));
  }
  const std::size_t iter = 10000;
  timings(dev, f, iter);
  timings(dev, g, iter);
}

void timings (DevicePtr_t dev, DifferentiableFunctionPtr_t f,
    const std::size_t iter)
{
  BasicConfigurationShooter cs (dev);
  std::cout << f->name() << '\n';
  const DifferentiableFunction& _f = *f;

  vector_t value = vector_t (f->outputSize ());
  clock_t value_elapsed = 0;
  boost::timer v_current;
  for (std::size_t i = 0; i < iter; i++) {
    ConfigurationPtr_t q = cs.shoot ();
    dev->currentConfiguration (*q);
    dev->computeForwardKinematics ();

    const clock_t begin_time = clock();
    _f (value, *q);
    value_elapsed += clock() - begin_time;
  }
  std::cout << "Value:\t" << value_elapsed << '\n';

  matrix_t jacobian = matrix_t (f->outputSize (), dev->numberDof ());
  clock_t jacobian_elapsed = 0;
  for (std::size_t i = 0; i < iter; i++) {
    ConfigurationPtr_t q = cs.shoot ();
    dev->currentConfiguration (*q);
    dev->computeForwardKinematics ();

    const clock_t begin_time = clock();
    _f.jacobian(jacobian, *q);
    jacobian_elapsed += clock() - begin_time;
  }
  std::cout << "Jacobian:\t" << jacobian_elapsed << '\n';
}

BOOST_AUTO_TEST_CASE (consistency) {
  DevicePtr_t device = createRobot ();
  JointPtr_t ee1 = device->getJointByName ("LLEG_5"),
             ee2 = device->getJointByName ("RLEG_5");
  Configuration_t goal = device->currentConfiguration ();
  BOOST_REQUIRE (device);
  BasicConfigurationShooter cs (device);

  /// The -1 factor is here because Orientation, Position, RelativeOrientation
  /// and RelativePosition express the relative position of Frame1 wrt Frame2, in
  /// Frame1. On the contrary, GenericTransformation express the relative transformation
  /// of Frame2 wrt Frame1, in Frame1.
  Transform3f Tid; Tid.setIdentity();
  check_consistent (device,
        Orientation::create ("Orientation"           , device, ee2, identity()),
        Orientation2::create("OrientationFromGeneric", device, ee2, Tid),
        -1);
        // Orientation::create (device, ee2, identity(), list_of(false)(true)(true))
  check_consistent (device,
        Position::create ("Position"           , device, ee2, vector3_t (0,0,0), vector3_t (0,0,0), identity ()),
        Position2::create("PositionFromGeneric", device, ee2, vector3_t (0,0,0), vector3_t (0,0,0)),
        -1);
        // Position::create (device, ee1, vector3_t (0,0,0), vector3_t (0,0,0), identity (), list_of(false)(true)(true))
  check_consistent (device,
        Position::create ("Position"           , device, ee2, vector3_t (1,0,0), vector3_t (0,0,0), identity ()),
        Position2::create("PositionFromGeneric", device, ee2, vector3_t (0,0,0), vector3_t (1,0,0)),
        -1);
        // Position::create (device, ee1, vector3_t (0,0,0), vector3_t (0,0,0), identity (), list_of(false)(true)(true))
  check_consistent (device,
        RelativeOrientation::create ("RelativeOrientation"           , device, ee1, ee2, identity ()),
        RelativeOrientation2::create("RelativeOrientationFromGeneric", device, ee1, ee2, Tid, Tid),
        -1);
        // RelativeOrientation::create (device, ee1, ee2, identity (), list_of(false)(true)(true))
  check_consistent (device,
        RelativePosition::create ("RelativePosition"           , device, ee1, ee2, vector3_t (0,0,0), vector3_t (0,0,0)),
        RelativePosition2::create("RelativePositionFromGeneric", device, ee1, ee2, vector3_t (0,0,0), vector3_t (0,0,0)),
        -1);
        // RelativePosition::create (device, ee1, ee2, vector3_t (0,0,0), vector3_t (0,0,0), list_of(false)(true)(true))

  ConfigurationPtr_t q2 = cs.shoot ();
  device->currentConfiguration (*q2);
  device->computeForwardKinematics ();
  Transform3f tf1 (device->getJointByName (device->name () + "_SO3")->
		   currentTransformation ());
  q2 = cs.shoot ();
  device->currentConfiguration (*q2);
  device->computeForwardKinematics ();
  Transform3f tf2 (device->getJointByName (device->name () + "_SO3")->
		   currentTransformation ());
  check_consistent (device,
      RelativeTransformation::create ("RelativeTransformation"           , device, ee1, ee2, tf1, tf2),
      RelativeTransformation2::create("RelativeTransformationFromGeneric", device, ee1, ee2, tf1, tf2));
}
