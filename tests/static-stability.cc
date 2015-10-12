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

#include "hpp/constraints/static-stability.hh"
#include "hpp/constraints/tools.hh"

#define BOOST_TEST_MODULE StaticStability
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
  const std::string& name = "";
  JointPtr_t joint, parent;
  std::string jointName = name + "_xyz";
  // Translation along xyz
  joint = objectFactory.createJointTranslation3 (transform3f_id());
  joint->name (jointName);
  for (size_type i = 0; i < 3; ++i) {
    joint->isBounded  (i, true);
    joint->lowerBound (i, -4);
    joint->upperBound (i, +4);
  }
  robot->rootJoint (joint);
  parent = joint;

  // joint SO3
  joint = objectFactory.createJointSO3 (transform3f_id());
  jointName = name + "_SO3";
  joint->name (jointName);
  parent->addChildJoint (joint);
  return joint;
}

JointPtr_t createRotationJoint (DevicePtr_t robot)
{
  const std::string& name = "";
  JointPtr_t joint;
  std::string jointName = name + "_rz";

  fcl::Transform3f mat; mat.setIdentity ();
  fcl::Matrix3f orient;
  orient (0,0) = 0; orient (0,1) = 0; orient (0,2) = 1;
  orient (1,0) = 0; orient (1,1) = 1; orient (1,2) = 0;
  orient (2,0) =-1; orient (2,1) = 0; orient (2,2) = 0;
  mat.setRotation (orient);

  // joint rz
  joint = objectFactory.createUnBoundedJointRotation (transform3f_id());
  joint->positionInParentFrame (mat);
  joint->name (jointName);
  robot->rootJoint (joint);
  return joint;
}

DevicePtr_t createRobot ()
{
  DevicePtr_t robot = Device::create ("StaticStability");
  JointPtr_t ff = createRotationJoint (robot);

  BodyPtr_t body;

  fcl::Transform3f mat; mat.setIdentity ();
  fcl::Matrix3f orient;
  orient (0,0) = 0; orient (0,1) = 0; orient (0,2) =-1;
  orient (1,0) = 0; orient (1,1) = 1; orient (1,2) = 0;
  orient (2,0) = 1; orient (2,1) = 0; orient (2,2) = 0;
  mat.setRotation (orient);

  // Translation along x
  JointPtr_t joint = objectFactory.createJointTranslation (transform3f_id());
  joint->name ("slider");
  ff->addChildJoint (joint);
  joint->positionInParentFrame (mat);
  joint->isBounded  (0, true);
  joint->lowerBound (0, -4);
  joint->upperBound (0, +4);
  body = objectFactory.createBody ();
  body->name ("slider");
  body->mass (1);
  joint->setLinkedBody (body);
  return robot;
}

StaticStabilityPtr_t createStaticStability (DevicePtr_t d, JointPtr_t j)
{
  StaticStability::Contacts_t cs;
  StaticStability::Contact_t c;

  c.joint1 = NULL; c.point1 = vector3_t (-1,0,0); c.normal1 = vector3_t (0,0,1);
  c.joint2 = j;    c.point2 = vector3_t (-1,0,0); c.normal2 = vector3_t (0,0,1);
  cs.push_back (c);

  c.joint1 = NULL; c.point1 = vector3_t (+1,0,0); c.normal1 = vector3_t (0,0,1);
  c.joint2 = j;    c.point2 = vector3_t (+1,0,0); c.normal2 = vector3_t (0,0,1);
  cs.push_back (c);

  CenterOfMassComputationPtr_t com = CenterOfMassComputation::create (d);
  com->add (d->rootJoint ());
  com->computeMass ();

  StaticStabilityPtr_t fptr = StaticStability::create ("test", d, cs, com);
  return fptr;
}

BOOST_AUTO_TEST_CASE (static_stability) {
  DevicePtr_t device = createRobot ();
  BOOST_REQUIRE (device);
  JointPtr_t rot = device->getJointByName ("_rz");
  JointPtr_t slider = device->getJointByName ("slider");

  Configuration_t c (3);
  c << 1, 0, 0;
  device->currentConfiguration (c);
  device->computeForwardKinematics ();
  // BOOST_MESSAGE ("\"rot\" initial transform:\n" << rot->currentTransformation ());
  // BOOST_MESSAGE ("\"rot\" current transform:\n" << rot->currentTransformation ());
  // BOOST_MESSAGE ("\"slider\" initial transform:\n" << slider->currentTransformation ());
  // BOOST_MESSAGE ("\"slider\" current transform:\n" << slider->currentTransformation ());

  BOOST_CHECK_MESSAGE (slider->currentTransformation ().isIdentity (),
      "This transform shoud be identity:\n" << slider->currentTransformation ());

  StaticStabilityPtr_t fptr  = createStaticStability (device, rot);
  StaticStabilityPtr_t fptrH = createStaticStabilityHard4 (device, rot);
  StaticStability& f  (*fptr);
  StaticStability& fH (*fptrH);
  const std::size_t nbC = 2;
  const std::size_t nbCH = 8;
  vector_t value (f.outputSize ());
  vector_t valueH (fH.outputSize ());
  matrix_t j (f.outputSize (), f.inputDerivativeSize ());
  std::list <Configuration_t> valid, invalid;

  valid.push_back ((Configuration_t (3) << 1,0,0).finished ());
  valid.push_back ((Configuration_t (3) << 1,0,0.5).finished ());
  valid.push_back ((Configuration_t (3) << 1,0,-0.5).finished ());
  valid.push_back ((Configuration_t (3) << 1,0,1).finished ());
  valid.push_back ((Configuration_t (3) << 1,0,-1).finished ());
  for (std::list<Configuration_t>::const_iterator it = valid.begin();
      it != valid.end (); ++it) {
    BOOST_MESSAGE ("Config " << it->transpose ());
    f (value, *it);
    BOOST_MESSAGE ("\"slider\" current transform:\n" << slider->currentTransformation ());
    BOOST_CHECK_MESSAGE (value.segment<6> (nbC).isZero (),
        "(I - phi * phi^+) * G =\n" << value.segment<6> (nbC).transpose());
    BOOST_CHECK_MESSAGE ((value.segment<nbCH> (0).array () > -Eigen::NumTraits<value_type>::dummy_precision()).all(),
        "Contact forces =\n" << value.segment<nbC> (0).transpose());
    fH (valueH, *it);
    BOOST_CHECK_MESSAGE (valueH.segment<6> (nbCH).isZero (),
        "(I - phi * phi^+) * G =\n" << valueH.segment<6> (nbCH).transpose());
    BOOST_CHECK_MESSAGE ((valueH.segment<nbCH> (0).array () > -Eigen::NumTraits<value_type>::dummy_precision()).all(),
        "Contact forces =\n" << valueH.segment<nbCH> (0).transpose());
  }

  BOOST_MESSAGE ("Starting invalid");
  invalid.push_back ((Configuration_t (3) << 0,1,0.5).finished ());
  invalid.push_back ((Configuration_t (3) << 1,0,2.5).finished ());
  invalid.push_back ((Configuration_t (3) << 1,0,-1.5).finished ());
  for (std::list<Configuration_t>::const_iterator it = invalid.begin();
      it != invalid.end (); ++it) {
    BOOST_MESSAGE ("Config " << it->transpose ());
    f (value, *it);
    BOOST_MESSAGE ("\"slider\" current transform:\n" << slider->currentTransformation ());
    BOOST_CHECK_MESSAGE (!value.segment<6> (nbC).isZero ()
        || (value.segment<nbC> (0).array () < - Eigen::NumTraits<value_type>::dummy_precision()).any(),
        "Should not be stable:\n"
        "(I - phi * phi^+) * G =\n" << value.segment<6> (nbC).transpose()
        << "\nContact forces =\n" << value.segment<nbC> (0).transpose()
        << "\nJoint transform:\n" << slider->currentTransformation ());

    fH (valueH, *it);
    BOOST_CHECK_MESSAGE (!valueH.segment<6> (nbCH).isZero ()
        || (valueH.segment<nbCH> (0).array () < - Eigen::NumTraits<value_type>::dummy_precision()).any(),
        "Should not be stable:\n"
        "(I - phi * phi^+) * G =\n" << valueH.segment<6> (nbCH).transpose()
        << "\nContact forces =\n" << valueH.segment<nbCH> (0).transpose()
        << "\nJoint transform:\n" << slider->currentTransformation ());
  }
}
