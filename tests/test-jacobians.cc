// Copyright (c) 2014, LAAS-CNRS
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
#include "hpp/constraints/relative-position.hh"
#include "hpp/constraints/relative-orientation.hh"
#include "hpp/constraints/relative-transformation.hh"
#include "hpp/constraints/static-stability.hh"

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
static fcl::Transform3f transform3f_id () { fcl::Transform3f T; T.setIdentity (); return T;}

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

StaticStabilityGravityPtr_t createStaticStability (DevicePtr_t d, JointPtr_t j)
{
  fcl::Vec3f x (1,0,0), y (0,1,0), z (0,0,1);
  fcl::Vec3f p[12];
  p[0] = fcl::Vec3f (-5,-5,0); p[1] = fcl::Vec3f (-5, 5,0); p[2] = fcl::Vec3f ( 5,-5,0);
  p[3] = fcl::Vec3f ( 5, 5,0); p[4] = fcl::Vec3f (-5, 5,0); p[6] = fcl::Vec3f ( 5,-5,0);
  p[6] = fcl::Vec3f ( 0, 0,1); p[7] = fcl::Vec3f (  1,0,1); p[8] = fcl::Vec3f (0,  1,1);
  p[9] = fcl::Vec3f ( 0, 0,0); p[10] = fcl::Vec3f (0.1,0,0); p[11] = fcl::Vec3f (0,0.1,0);
  fcl::TriangleP f1 (p[0],p[1],p[2]),
                 f2 (p[3],p[4],p[5]),
                 th (p[6],p[7],p[8]),
                 o  (p[9],p[10],p[11]);
  StaticStabilityGravityPtr_t fptr = StaticStabilityGravity::create (d, j);
  StaticStabilityGravity& f = *fptr;
  f.addObjectTriangle (o);
  f.addFloorTriangle (th);
  f.addFloorTriangle (f1);
  f.addFloorTriangle (f2);
  return fptr;
}

BOOST_AUTO_TEST_CASE (triangle) {
  /// First test Triangle class
  fcl::Vec3f x (1,0,0), y (0,1,0), z (0,0,1);
  fcl::Vec3f p[9];
  p[0] = fcl::Vec3f (0,0,0);
  p[1] = fcl::Vec3f (1,0,0);
  p[2] = fcl::Vec3f (0,1,0);
  p[3] = fcl::Vec3f (1,1,1);
  p[4] = fcl::Vec3f (0.2,0.2,0);
  p[5] = fcl::Vec3f (1,1,0);
  p[6] = fcl::Vec3f (-1,-1,1);
  Triangle t (p[0],p[1],p[2]);
  BOOST_CHECK_MESSAGE ((t.normal () - z).isZero (), "Norm of triangle is wrong");
  BOOST_CHECK_MESSAGE ((t.center () - (x+y)/3).isZero (), "Center of triangle is wrong");
  BOOST_CHECK_MESSAGE (std::abs (t.planeXaxis ().dot (z)) < 1e-8, "X axis of triangle is wrong");
  BOOST_CHECK_MESSAGE (std::abs (t.planeYaxis ().dot (z)) < 1e-8, "Y axis of triangle is wrong");
  BOOST_CHECK_MESSAGE (std::abs (t.planeYaxis ().dot (t.planeXaxis ())) < 1e-8, "X axis of triangle is wrong");
  BOOST_CHECK_MESSAGE ((t.intersection (p[6], y-z) + x).isZero (), "Wrong intersection of triangle and line is wrong");
  BOOST_CHECK_MESSAGE (t.isInside (p[4]), "This point is inside");
  BOOST_CHECK_MESSAGE (!t.isInside (p[5]), "This point is outside");
  BOOST_CHECK_MESSAGE (std::abs (t.distance (p[0])) < 1e-8, "Distance to triangle is wrong");
  BOOST_CHECK_MESSAGE (std::abs (t.distance (p[1])) < 1e-8, "Distance to triangle is wrong");
  BOOST_CHECK_MESSAGE (std::abs (t.distance (p[2])) < 1e-8, "Distance to triangle is wrong");
  BOOST_CHECK_MESSAGE (std::abs (t.distance (p[4]) + 0.2) < 1e-8, "Distance to triangle is wrong");
  BOOST_CHECK_MESSAGE (std::abs (t.distance (p[5]) - 0.5 * std::sqrt(2)) < 1e-8, "Distance to triangle is wrong");
  BOOST_CHECK_MESSAGE (std::abs (t.distance (t.intersection (p[6], z)) - std::sqrt(2)) < 1e-8, "Distance to triangle is wrong");
}

BOOST_AUTO_TEST_CASE (static_stability) {
  DevicePtr_t device = createRobot ();
  JointPtr_t ee1 = device->getJointByName ("LLEG_5");
  BOOST_REQUIRE (device);

  StaticStabilityGravityPtr_t fptr = createStaticStability (device, ee1);
  StaticStabilityGravity& f = *fptr;
  vector_t value (f.outputSize ());
  f (value, device->currentConfiguration ());
  matrix_t j (f.outputSize (), f.inputDerivativeSize ());
  f.jacobian (j, device->currentConfiguration ());
}

BOOST_AUTO_TEST_CASE (jacobian) {
  DevicePtr_t device = createRobot ();
  JointPtr_t ee1 = device->getJointByName ("LLEG_5"),
             ee2 = device->getJointByName ("RLEG_5");
  BOOST_REQUIRE (device);
  BasicConfigurationShooter cs (device);

  /// Create the constraints
  typedef DifferentiableFunction DF;
  typedef std::pair <std::string, DifferentiableFunctionPtr_t> DFptr;
  typedef std::list <DFptr> DFs;
  DFs functions;
  functions.push_back ( DFptr (
        "Orientation",
        Orientation::create (device, ee2, identity())
      ));
  functions.push_back ( DFptr (
        "Orientation with mask (0,1,1)",
        Orientation::create (device, ee2, identity(), list_of(false)(true)(true))
      ));
  functions.push_back ( DFptr (
        "Position",
        Position::create (device, ee1, vector3_t (0,0,0),
          vector3_t (0,0,0), identity ())
      ));
  functions.push_back ( DFptr (
        "Position with mask (0,1,1)",
        Position::create (device, ee1, vector3_t (0,0,0),
          vector3_t (0,0,0), identity (), list_of(false)(true)(true))
      ));
  functions.push_back ( DFptr (
        "RelativeOrientation",
        RelativeOrientation::create (device, ee1, ee2, identity ())
      ));
  functions.push_back ( DFptr (
        "RelativeOrientation with mask (0,1,1)",
        RelativeOrientation::create (device, ee1, ee2, identity (), list_of(false)(true)(true))
      ));
  functions.push_back ( DFptr (
        "RelativePosition",
        RelativePosition::create (device, ee1, ee2, vector3_t (0,0,0),
          vector3_t (0,0,0))
      ));
  functions.push_back ( DFptr (
        "RelativePosition with mask (0,1,1)",
        RelativePosition::create (device, ee1, ee2, vector3_t (0,0,0),
          vector3_t (0,0,0), list_of(false)(true)(true))
      ));
  functions.push_back ( DFptr (
        "RelativeTransformation",
        RelativeTransformation::create (device, ee1, ee2, transform3f_id ())
      ));
  functions.push_back ( DFptr (
        "StaticStabilityGravity",
        createStaticStability (device, ee1)
      ));

  ConfigurationPtr_t q1, q2 = cs.shoot ();
  vector_t value1, value2, dvalue, error;
  vector_t errorNorm (MAX_NB_ERROR);
  vector_t dq (device->numberDof ()); dq.setZero ();
  matrix_t jacobian;
  for (DFs::iterator fit = functions.begin(); fit != functions.end(); ++fit) {
    DF& f = *(fit->second);
    value1 = vector_t (f.outputSize ());
    value2 = vector_t (f.outputSize ());
    errorNorm.setZero ();
    jacobian = matrix_t (f.outputSize (), device->numberDof ());
    for (size_t i = 0; i < NUMBER_JACOBIAN_CALCULUS; i++) {
      q1 = cs.shoot ();
      f (value1, *q1);
      jacobian.setZero ();
      f.jacobian (jacobian, *q1);
      // We check the jacobian for each DOF.
      for (int idof = 0; idof < device->numberDof (); idof++){
        dvalue = jacobian.col (idof);
        // dq = (0,...,0,1,0,...,0), the 1 being at the rank idof.
        // Check that ( e(q1 + eps*dq) - e(q1) / eps) -> jacobian * dq
        size_t i_error;
        dq[idof] = 10 * DQ_MAX;
        for (i_error = 0; i_error < MAX_NB_ERROR; i_error++) {
          //dq[idof] = DQ_MAX * std::pow (10, - i_error);
          dq[idof] = dq[idof] / 10;
          hpp::model::integrate (device, *q1, dq, *q2);
          f (value2, *q2);
          error = value2 - value1 - dq[idof] * dvalue;
          errorNorm [i_error] = error.norm ();
          if (errorNorm [i_error] < 0.5 * dq[idof] * dq[idof] * HESSIAN_MAXIMUM_COEF)
            break;
        }
        BOOST_CHECK_MESSAGE (i_error != MAX_NB_ERROR,
              "Constraint " << fit->first << ": error norm " << errorNorm [MAX_NB_ERROR - 1] / dq[idof]
              << ", dof " << idof << ", config " << i << ".");
        dq(idof) = 0;
      }
    }
  }
}
