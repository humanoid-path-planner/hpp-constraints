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
#include <hpp/core/basic-configuration-shooter.hh>

#include "hpp/constraints/position.hh"
#include "hpp/constraints/orientation.hh"
#include "hpp/constraints/relative-position.hh"
#include "hpp/constraints/relative-orientation.hh"

#define BOOST_TEST_MODULE hpp_constraints
#include <boost/test/included/unit_test.hpp>

#include <stdlib.h>
#include <math.h>

using hpp::model::Device;
using hpp::model::DevicePtr_t;
using hpp::model::JointPtr_t;
using hpp::model::JointSO3;
using hpp::model::JointAnchor;
using hpp::model::JointRotation;
using hpp::model::JointTranslation;

using hpp::core::ConfigurationPtr_t;
using hpp::core::BasicConfigurationShooter;

using boost::assign::list_of;

using namespace hpp::constraints;

const static size_t NUMBER_JACOBIAN_CALCULUS = 5;
const static double HESSIAN_MAXIMUM_COEF = 1e2;
const static double DQ_MAX = 1e-2;
const static size_t MAX_NB_ERROR = 5;

static matrix3_t identity () { matrix3_t R; R.setIdentity (); return R;}

DevicePtr_t createDevice (JointPtr_t& endEffector1, JointPtr_t& endEffector2) {
  DevicePtr_t robot = Device::create("test");
  JointPtr_t xJoint = new JointTranslation(fcl::Transform3f());
  xJoint->isBounded(0,1);
  xJoint->lowerBound(0,-3.);
  xJoint->upperBound(0,3.);
  JointPtr_t yJoint = new JointTranslation
    (fcl::Transform3f(fcl::Quaternion3f (sqrt (2)/2, 0, 0, sqrt(2)/2)));
  yJoint->isBounded(0,1);
  yJoint->lowerBound(0,-3.);
  yJoint->upperBound(0,3.);
  JointPtr_t zJoint = new JointTranslation
    (fcl::Transform3f(fcl::Quaternion3f (sqrt (2)/2, 0, sqrt(2)/2, 0)));
  zJoint->isBounded(0,1);
  zJoint->lowerBound(0,-3.);
  zJoint->upperBound(0,3.);
  JointPtr_t so3Joint = new JointSO3(fcl::Transform3f());

  robot->rootJoint (xJoint);
  xJoint->addChildJoint (yJoint);
  xJoint->addChildJoint (zJoint);
  zJoint->addChildJoint (so3Joint);

  /// Adding two arms
  JointPtr_t arms[2];
  JointPtr_t hands[2];
  JointPtr_t endEffectors[2];
  fcl::Matrix3f xtoz; xtoz.setEulerZYX (0, -M_PI/2, 0);
  for (size_t i = 0; i < 2; i++) {
    arms[i] = new JointRotation (fcl::Transform3f (
          xtoz,
          fcl::Vec3f ( 2*i - 1, 0, 0)
          ));
    so3Joint->addChildJoint (arms[i]);
    hands[i] = new JointRotation (fcl::Transform3f (
          xtoz,
          fcl::Vec3f ( 1, 0, 0)
          ));
    arms[i]->addChildJoint (hands[i]);
    endEffectors[i] = new JointAnchor (
          fcl::Vec3f ( 1, 0, 0)
        );
    hands[i]->addChildJoint (endEffectors [i]);
  }
  endEffector1 = endEffectors[0];
  endEffector2 = endEffectors[1];

  return robot;
}

void shootDQ (vector_t& dq) {
  for (int i = 0; i < dq.size (); i++) {
    dq [i] = (- 1 + 2 * ((double) rand () / RAND_MAX)) * DQ_MAX;
  }
}

BOOST_AUTO_TEST_CASE (jacobian) {
  JointPtr_t ee1, ee2;
  DevicePtr_t device = createDevice (ee1, ee2);
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

  ConfigurationPtr_t q1, q2 = cs.shoot ();
  vector_t value1, value2, dvalue, error;
  vector_t errorNorm (MAX_NB_ERROR);
  vector_t dq (device->numberDof ()); dq.setZero ();
  matrix_t jacobian;
  for (DFs::iterator fit = functions.begin(); fit != functions.end(); fit++) {
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
              "Constraint " << fit->first << ": error norm " << errorNorm [MAX_NB_ERROR - 1]
              << ", dof " << idof << ", config " << i << ".");
        dq(idof) = 0;
      }
    }
  }
}
