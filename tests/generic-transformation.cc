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

using hpp::model::Configuration_t;
using hpp::model::ConfigurationPtr_t;
using hpp::model::Device;
using hpp::model::DevicePtr_t;
using hpp::model::JointPtr_t;
using hpp::model::JointVector_t;
using hpp::model::BodyPtr_t;
using hpp::model::Transform3f;

using namespace hpp::constraints;

const static size_t NUMBER_JACOBIAN_CALCULUS = 10;

static matrix3_t identity () { matrix3_t R; R.setIdentity (); return R;}

hpp::model::ObjectFactory objectFactory;

matrix3_t exponential (const vector3_t& aa)
{
  matrix3_t R, xCross;
  xCross.setZero();
  xCross(1, 0) = + aa(2); xCross(0, 1) = - aa(2);
  xCross(2, 0) = - aa(1); xCross(0, 2) = + aa(1);
  xCross(2, 1) = + aa(0); xCross(1, 2) = - aa(0);
  R.setIdentity();
  value_type theta = aa.norm();
  if (theta < 1e-6) {
    R += xCross;
    R += 0.5 * xCross.transpose() * xCross;
  } else {
    R += sin(theta) / theta * xCross;
    R += 2 * std::pow(sin(theta/2),2) / std::pow(theta,2) * xCross.transpose() * xCross;
  }
  return R;
}

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
void timings2 (DevicePtr_t dev,
    const std::vector<ConfigurationPtr_t>& configs,
    DifferentiableFunctionPtr_t f);

struct ProportionalCompare {
  const value_type alpha;
  ProportionalCompare (value_type _alpha = 1) : alpha (_alpha) {}

  void value(const vector_t& value1, const vector_t& value2) const {
    vector_t d = value2 - alpha * value1;
    //std::cout << "----------" << std::endl;
    //std::cout << value1.transpose() << std::endl;
    //std::cout << value2.transpose() << std::endl;
    //std::cout << "--" << std::endl;
    //std::cout << exponential(value1.tail<3>()) * R << std::endl;
    //std::cout << exponential(value2.tail<3>()) * transpose(R) << std::endl;
    //std::cout << "----------" << std::endl;
    //std::cout << d.transpose() << std::endl;
    BOOST_CHECK_MESSAGE(value1.isApprox(alpha*value2, 1e-10),
			"Value not matching. Norm of value1 is "
			<< value1.norm () <<  ", norm of value2 is "
			<< value2.norm () << ", norm of error is "
			<< d.norm());
  }

  void jacobian(const matrix_t& jacobian1, const matrix_t& jacobian2) const {
    matrix_t diffJ = jacobian2 - alpha*jacobian1;
    // std::cout << diffJ.norm() << std::endl;
    BOOST_CHECK_MESSAGE(jacobian1.isApprox(alpha*jacobian2, 1e-10),
			"Jacobian not matching. Norm of J1 is "
			<< jacobian1.norm () << ", norm of J2 is "
			<< jacobian2.norm () << ", norm of error is "
			<< diffJ.norm());
  }
};

struct TransformationCompare {
  const matrix3_t R;
  const size_type rowTr, sizeTr;
  TransformationCompare (const matrix3_t _R, size_type _rowTr, size_type _sizeTr)
    : R (_R), rowTr(_rowTr), sizeTr(_sizeTr) {}

  void value(const vector_t& value1, const vector_t& value2) const {
    vector_t dtr  = R.middleRows(rowTr, sizeTr) * value2.segment(rowTr,sizeTr) + value1.segment(rowTr,sizeTr);
    BOOST_CHECK_MESSAGE(dtr.isZero(1e-10),
			"Value not matching. Norm of value1 is "
			<< value1.segment(rowTr,sizeTr).norm () <<  ", norm of value2 is "
			<< value2.segment(rowTr,sizeTr).norm () << ", norm of error is "
			<< dtr.norm());
    //matrix3_t A = exponential(value1.tail<3>()) * transpose(R);
    //matrix3_t B = exponential(value2.tail<3>()) * R;
    matrix3_t A = exponential(value1.tail<3>());
    matrix3_t B = exponential(value2.tail<3>());
    matrix_t diff = A.transpose() - B;
    // std::cout << diff << std::endl;
    BOOST_CHECK_MESSAGE(diff.isZero(1e-10),
			"Value not matching. Norm of error is "
			<< diff.norm());
  }

  void jacobian(const matrix_t& jacobian1, const matrix_t& jacobian2) const {
    matrix_t diffJtr = R.middleRows(rowTr, sizeTr) * jacobian2.middleRows(rowTr, sizeTr) + jacobian1.middleRows(rowTr, sizeTr);
    //matrix_t diffJrot = jacobian2 - jacobian1;
    // std::cout << diffJ.norm() << std::endl;
    BOOST_CHECK_MESSAGE(diffJtr.isZero(1e-10),
			"Jacobian not matching. Norm of J1 is "
			<< jacobian1.middleRows(rowTr, sizeTr).norm () << ", norm of J2 is "
			<< jacobian2.middleRows(rowTr, sizeTr).norm () << ", norm of error is "
			<< diffJtr.norm());
    std::cout << "Rotation part of jacobian not checked." << std::endl;
  }
};

template <typename Compare>
void check_consistent (DevicePtr_t dev,
    DifferentiableFunctionPtr_t f, DifferentiableFunctionPtr_t g,
    const Compare comp = Compare())
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
    comp.value (value1, value2);
    f->jacobian (jacobian1, *q);
    g->jacobian (jacobian2, *q);
    comp.jacobian (jacobian1, jacobian2);
  }
  const std::size_t iter = 10000;
  // timings(dev, f, iter);
  // timings(dev, g, iter);
  //std::vector<ConfigurationPtr_t> cfgs(iter);
  //for (size_t i = 0; i < iter; i++) cfgs[i] = cs.shoot();
  //timings2(dev, cfgs, f);
  //timings2(dev, cfgs, g);
}

void timings (DevicePtr_t dev, DifferentiableFunctionPtr_t f,
    const std::size_t iter)
{
  BasicConfigurationShooter cs (dev);
  std::cout << "=======================\n" << f->name() << '\n';
  const DifferentiableFunction& _f = *f;

  vector_t value = vector_t (f->outputSize ());
  clock_t value_elapsed = 0;
  for (std::size_t i = 0; i < iter; i++) {
    ConfigurationPtr_t q = cs.shoot ();
    dev->currentConfiguration (*q);
    dev->computeForwardKinematics ();

    const clock_t begin_time = clock();
    _f (value, *q);
    value_elapsed += clock() - begin_time;
  }
  std::cout << "Value   :\t" << value_elapsed << '\n';

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

void timings2 (DevicePtr_t dev,
    const std::vector<ConfigurationPtr_t>& configs,
    DifferentiableFunctionPtr_t f)
{
  std::cout << "=======================\n" << f->name() << '\n';
  const DifferentiableFunction& _f = *f;
  const std::size_t iter = configs.size();

  vector_t value = vector_t (f->outputSize ());
  matrix_t jacobian = matrix_t (f->outputDerivativeSize (), f->inputDerivativeSize ());
  
  clock_t elapsed = 0;
  for (std::size_t i = 0; i < iter; i++) {
    dev->currentConfiguration (*configs[i]);
    dev->computeForwardKinematics ();
    const clock_t begin_time = clock();
    _f (value, *configs[i]);
    _f.jacobian(jacobian, *configs[i]);
    elapsed += clock() - begin_time;
  }
  std::cout << "Time   :\t" << elapsed << '\n';
}

BOOST_AUTO_TEST_CASE (consistency) {
  DevicePtr_t device = createRobot ();
  JointPtr_t ee1 = device->getJointByName ("LLEG_5"),
             ee2 = device->getJointByName ("RLEG_5");
  Configuration_t goal = device->currentConfiguration ();
  BOOST_REQUIRE (device);
  BasicConfigurationShooter cs (device);
  std::vector<bool> mask011 (3, true); mask011[0] = false;

  device->currentConfiguration (*cs.shoot ());
  device->computeForwardKinematics ();
  Transform3f tf1 (ee1->currentTransformation ());
  Transform3f tf2 (ee2->currentTransformation ());

  /// The -1 factor is here because Orientation, Position, RelativeOrientation
  /// and RelativePosition express the relative position of Frame1 wrt Frame2, in
  /// Frame1. On the contrary, GenericTransformation express the relative transformation
  /// of Frame2 wrt Frame1, in Frame1.
  Transform3f Tid; Tid.setIdentity();
  check_consistent (device,
        deprecated::Orientation::create ("Orientation"           , device, ee2, tf2.getRotation()),
        Orientation::create             ("OrientationFromGeneric", device, ee2, tf2.getRotation()),
        ProportionalCompare(-1));
  check_consistent (device,
        deprecated::Orientation::create ("Orientation"           , device, ee2, tf2.getRotation(), mask011),
        Orientation::create             ("OrientationFromGeneric", device, ee2, tf2.getRotation(), mask011),
        ProportionalCompare(-1));
  check_consistent (device,
        deprecated::Position::create ("Position"           , device, ee2, tf2.getTranslation(), tf1.getTranslation(), transpose(tf1.getRotation())),
        Position::create             ("PositionFromGeneric", device, ee2, tf2, tf1),
        ProportionalCompare(-1));
  check_consistent (device,
        deprecated::Position::create ("Position"           , device, ee2, tf2.getTranslation(), tf1.getTranslation(), transpose(tf1.getRotation()), mask011),
        Position::create             ("PositionFromGeneric", device, ee2, tf2, tf1, mask011),
        ProportionalCompare(-1));
  check_consistent (device,
        deprecated::RelativeOrientation::create ("RelativeOrientation"           , device, ee1, ee2, tf1.getRotation()),
        RelativeOrientation::create             ("RelativeOrientationFromGeneric", device, ee1, ee2, tf1.getRotation()),
        ProportionalCompare(-1));
  check_consistent (device,
        deprecated::RelativeOrientation::create ("RelativeOrientation"           , device, ee1, ee2, tf1.getRotation(), mask011),
        RelativeOrientation::create             ("RelativeOrientationFromGeneric", device, ee1, ee2, tf1.getRotation(), mask011),
        ProportionalCompare(-1));
  check_consistent (device,
        deprecated::RelativePosition::create ("RelativePosition"           , device, ee1, ee2, tf1.getTranslation(), tf2.getTranslation()),
        RelativePosition::create             ("RelativePositionFromGeneric", device, ee1, ee2, tf1.getTranslation(), tf2.getTranslation()),
        ProportionalCompare(-1));

  std::vector < bool > mask (6, true); mask [0] = mask [1] = false;
  check_consistent (device,
        deprecated::Orientation::create ("Orientation"           , device, ee2, tf2.getRotation()),
        Orientation::create             ("OrientationFromGeneric", device, ee2, tf2.getRotation()),
        ProportionalCompare(-1));
  check_consistent (device,
        deprecated::Position::create ("Position"           , device, ee2, tf2.getTranslation(), tf1.getTranslation(), transpose(tf1.getRotation())),
        Position::create             ("PositionFromGeneric", device, ee2, tf2.getTranslation(), tf1),
        ProportionalCompare(-1));
  check_consistent (device, deprecated::RelativeTransformation::create
		    ("RelativeTransformation", device, ee1, ee2, Tid, Tid, mask),
		    RelativeTransformation::create
		    ("RelativeTransformationFromGeneric", device, ee1, ee2, Tid, Tid, mask),
        ProportionalCompare(1));
  check_consistent (device,
      deprecated::RelativeTransformation::create ("RelativeTransformation"           , device, ee1, ee2, tf1, tf2),
      RelativeTransformation::create             ("RelativeTransformationFromGeneric", device, ee1, ee2, tf1, tf2),
      ProportionalCompare(1));
  check_consistent (device,
		    deprecated::Transformation::create ("Transformation", device, ee1, identity()),
		    Transformation::create ("TransformationFromGeneric", device, ee1, identity()),
        ProportionalCompare(-1));

  fcl::Matrix3f R = transpose(tf2.getRotation());
  fcl::Transform3f tf2inv = inverse(tf2);
  check_consistent (device,
		    deprecated::Transformation::create ("Transformation", device, ee1, tf2),
		    Transformation::create ("TransformationFromGeneric", device, ee1, tf2),
        TransformationCompare (tf2.getRotation(), 0, 3));
  //check_consistent (device,
		    //deprecated::Transformation::create ("Transformation", device, ee1, tf2, mask),
		    //Transformation::create ("TransformationFromGeneric", device, ee1, tf2, mask));
  check_consistent (device,
        deprecated::RelativeOrientation::create ("RelativeOrientation"           , device, ee1, ee2, tf1.getRotation()),
        RelativeOrientation::create             ("RelativeOrientationFromGeneric", device, ee1, ee2, tf1.getRotation()),
        ProportionalCompare(-1));
}
