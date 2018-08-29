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

#include <pinocchio/parsers/sample-models.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/simple-device.hh>

#include "hpp/constraints/generic-transformation.hh"
#include "hpp/constraints/symbolic-function.hh"
#include "hpp/constraints/convex-shape-contact.hh"
#include "hpp/constraints/static-stability.hh"
#include "hpp/constraints/configuration-constraint.hh"
#include "hpp/constraints/differentiable-function-stack.hh"
#include "hpp/constraints/tools.hh"

// The following 6 includes are deprecated and will be removed
// #include "hpp/constraints/position.hh"
// #include "hpp/constraints/orientation.hh"
// #include "hpp/constraints/transformation.hh"
// #include "hpp/constraints/relative-position.hh"
// #include "hpp/constraints/relative-orientation.hh"
// #include "hpp/constraints/relative-transformation.hh"

#define BOOST_TEST_MODULE hpp_constraints
#include <boost/test/included/unit_test.hpp>

#include <stdlib.h>
#include <limits>
#include <math.h>

using hpp::pinocchio::Configuration_t;
using hpp::pinocchio::ConfigurationPtr_t;
using hpp::pinocchio::Device;
using hpp::pinocchio::DevicePtr_t;
using hpp::pinocchio::JointPtr_t;
using hpp::pinocchio::JointVector_t;
using hpp::pinocchio::BodyPtr_t;

using std::numeric_limits;
using boost::assign::list_of;

typedef std::vector<bool> BoolVector_t;

using namespace hpp::constraints;

const static size_t NUMBER_JACOBIAN_CALCULUS = 5;
const static double HESSIAN_MAXIMUM_COEF = 1e1;
const static double DQ_MAX = 1e-2;
const static size_t MAX_NB_ERROR = 5;

static matrix3_t I3 = matrix3_t::Identity();
static vector3_t zero = vector3_t::Identity();
static se3::SE3 MId = se3::SE3::Identity();
se3::SE3 toSE3(const matrix3_t& R) { return se3::SE3(R, zero); }
se3::SE3 toSE3(const vector3_t& t) { return se3::SE3(I3, t); }

class BasicConfigurationShooter
{
public:
  BasicConfigurationShooter (const DevicePtr_t& robot) : robot_ (robot)
  {
  }
  virtual ConfigurationPtr_t shoot () const
  {
    ConfigurationPtr_t config (new Configuration_t (robot_->configSize ()));
    *config = se3::randomConfiguration(robot_->model());
    return config;
  }
private:
  const DevicePtr_t& robot_;
}; // class BasicConfigurationShooter

DevicePtr_t createRobot ()
{
  DevicePtr_t robot =
    hpp::pinocchio::humanoidSimple("test", true,
        (Device::Computation_t) (Device::JOINT_POSITION | Device::JACOBIAN));
  robot->rootJoint()->lowerBound(0, -1);
  robot->rootJoint()->lowerBound(1, -1);
  robot->rootJoint()->lowerBound(2, -1);
  robot->rootJoint()->upperBound(0,  1);
  robot->rootJoint()->upperBound(1,  1);
  robot->rootJoint()->upperBound(2,  1);
  return robot;
}

ConvexShapeContactPtr_t createConvexShapeContact_punctual (DevicePtr_t d, JointPtr_t j, std::string name)
{
  /** Floor = penta + square
   *     +
   *    / \
   *   /   +-----+
   *  +    |     |
   *   \   +-----+
   *    \ /
   *     +
   *
   *  Object = point
  **/
  std::vector <vector3_t> square(4);
  square[0] = vector3_t ( 5, 5,0); square[1] = vector3_t ( 5,-5,0);
  square[2] = vector3_t ( 0,-5,0); square[3] = vector3_t ( 0, 5,0);
  ConvexShape sCs (square);

  std::vector <vector3_t> penta(5);
  penta[0] = vector3_t ( 0, 5,0); penta[1] = vector3_t ( 0,-5,0);
  penta[2] = vector3_t (-2,-6,0); penta[3] = vector3_t (-5, 0,0);
  penta[4] = vector3_t (-2, 6,0);
  ConvexShape pCs (penta);

  std::vector <vector3_t> point(1, vector3_t (0,0,0));
  ConvexShape pointCs (point, j);

  ConvexShapeContactPtr_t fptr = ConvexShapeContact::create (name, d);
  ConvexShapeContact& f = *fptr;
  f.addObject (pointCs);
  f.addFloor (sCs);
  f.addFloor (pCs);
  return fptr;
}

ConvexShapeContactPtr_t createConvexShapeContact_convex (DevicePtr_t d, JointPtr_t j, std::string name)
{
  /** Floor = penta + square
   *     +
   *    / \
   *   /   +-----+
   *  +    |     |
   *   \   +-----+
   *    \ /
   *     +
   *
   *  Object = point
   *  +--+
   *  |   \
   *  +----+
  **/
  std::vector <vector3_t> square(4);
  square[0] = vector3_t ( 5, 5,0); square[1] = vector3_t ( 5,-5,0);
  square[2] = vector3_t ( 0,-5,0); square[3] = vector3_t ( 0, 5,0);
  ConvexShape sCs (square);

  std::vector <vector3_t> penta(5);
  penta[0] = vector3_t ( 0, 5,0); penta[1] = vector3_t ( 0,-5,0);
  penta[2] = vector3_t (-2,-6,0); penta[3] = vector3_t (-5, 0,0);
  penta[4] = vector3_t (-2, 6,0);
  ConvexShape pCs (penta);

  std::vector <vector3_t> trapeze(4);
  trapeze[0] = vector3_t (-0.1, 0.1,0); trapeze[1] = vector3_t ( 0.1, 0.1,0);
  trapeze[2] = vector3_t ( 0.2,-0.1,0); trapeze[3] = vector3_t (-0.1,-0.1,0);
  ConvexShape tCs (trapeze, j);

  ConvexShapeContactPtr_t fptr = ConvexShapeContact::create (name, d);
  ConvexShapeContact& f = *fptr;
  f.addObject (tCs);
  f.addFloor (sCs);
  f.addFloor (pCs);
  return fptr;
}

ConvexShapeContactPtr_t createConvexShapeContact_triangles (DevicePtr_t d, JointPtr_t j, std::string name)
{
  vector3_t x (1,0,0), y (0,1,0), z (0,0,1);
  vector3_t p[12];
  p[0] = vector3_t (-5,-5,0); p[1] = vector3_t (-5, 5,0); p[2] = vector3_t ( 5,-5,0);
  p[3] = vector3_t ( 5, 5,0); p[4] = vector3_t (-5, 5,0); p[5] = vector3_t ( 5,-5,0);
  p[6] = vector3_t ( 0, 0,1); p[7] = vector3_t (  1,0,1); p[8] = vector3_t (0,  1,1);
  p[9] = vector3_t ( 0, 0,0); p[10] = vector3_t (0.1,0,0); p[11] = vector3_t (0,0.1,0);
  std::vector<vector3_t> f1 = list_of (p[0])(p[1])(p[2]),
                         f2 = list_of (p[3])(p[4])(p[5]),
                         th = list_of (p[6])(p[7])(p[8]),
                         o  = list_of (p[9])(p[10])(p[11]);
  ConvexShapeContactPtr_t fptr = ConvexShapeContact::create (name, d);
  ConvexShapeContact& f = *fptr;
  f.addObject (ConvexShape (o, j));
  f.addFloor (ConvexShape (th, JointPtr_t()));
  f.addFloor (ConvexShape (f1, JointPtr_t()));
  f.addFloor (ConvexShape (f2, JointPtr_t()));
  return fptr;
}

BOOST_AUTO_TEST_CASE (triangle) {
  /// First test ConvexShape class (as a triangle)
  vector3_t x (1,0,0), y (0,1,0), z (0,0,1);
  vector3_t p[9];
  p[0] = vector3_t (0,0,0);
  p[1] = vector3_t (1,0,0);
  p[2] = vector3_t (0,1,0);
  p[3] = vector3_t (1,1,1);
  p[4] = vector3_t (0.2,0.2,0);
  p[5] = vector3_t (1,1,0);
  p[6] = vector3_t (-1,-1,1);
  ConvexShape t (p[0],p[1],p[2]);
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

template<bool forward>
void checkJacobianDiffIsZero(const std::string& name, const matrix_t& diff, const value_type& eps)
{
  size_type row, col;
  value_type maxError = diff.cwiseAbs().maxCoeff(&row,&col);
  BOOST_CHECK_MESSAGE(maxError < /* HESSIAN_MAXIMUM_COEF * */ eps,
      "Jacobian of " << name << " seems wrong (" << (forward? "forward":"central") <<  "). DOF " << col << " at row " << row << ": "
      << maxError << " > " << /* HESSIAN_MAXIMUM_COEF << " * " << */ eps
      );
}

BOOST_AUTO_TEST_CASE (jacobian) {
  DevicePtr_t device = createRobot ();
  JointPtr_t ee1 = device->getJointByName ("lleg5_joint"),
             ee2 = device->getJointByName ("rleg5_joint");
  Configuration_t goal = device->currentConfiguration ();
  BOOST_REQUIRE (device);
  BasicConfigurationShooter cs (device);

  device->currentConfiguration (*cs.shoot ());
  device->computeForwardKinematics ();
  Transform3f tf1 (ee1->currentTransformation ());
  Transform3f tf2 (ee2->currentTransformation ());

  /// Create the constraints
  typedef std::list <DifferentiableFunctionPtr_t> DFs;
  std::vector<bool> mask011 (3, true); mask011[0] = false;
  DFs functions;
  functions.push_back (
      ConfigurationConstraint::create ("Configuration", device, goal));
  // /*
  functions.push_back (
      Orientation::create ("Orientation", device, ee2, toSE3(tf2.rotation())));
  functions.push_back (
      Orientation::create ("Orientation with mask (0,1,1)", device, ee1, toSE3(tf1.rotation()), mask011));
  functions.push_back (
      Position::create ("Position", device, ee1, toSE3(tf1.translation()), tf2));
  functions.push_back (
      Position::create ("Position with mask (0,1,1)", device, ee1, toSE3(tf1.translation()), tf2, mask011));
  functions.push_back (
      RelativeOrientation::create ("RelativeOrientation", device, ee1, ee2, tf1, tf2));
  functions.push_back (
      RelativeOrientation::create ("RelativeOrientation with mask (0,1,1)", device, ee1, ee2, tf1, tf2, mask011));
  functions.push_back (
      RelativePosition::create ("RelativePosition", device, ee1, ee2, tf1, tf2));
  functions.push_back (
      RelativePosition::create ("RelativePosition with mask (0,1,1)", device, ee1, ee2, tf1, tf2, mask011));

  device->currentConfiguration (*cs.shoot ());
  device->computeForwardKinematics ();
  tf1 = ee1->currentTransformation ();
  tf2 = ee2->currentTransformation ();

  functions.push_back (
      Transformation::create ("Transformation", device, ee1, tf1));
  functions.push_back (
      RelativeTransformation::create ("RelativeTransformation", device, ee1, ee2, tf1, tf2));
  functions.push_back (
      createConvexShapeContact_triangles (device, ee1, "ConvexShapeContact triangle"));
  functions.push_back (
      createConvexShapeContact_punctual (device, ee1, "ConvexShapeContact punctual"));
  functions.push_back (
      createConvexShapeContact_convex (device, ee1, "ConvexShapeContact convex"));

  // DifferentiableFunctionSet
  DifferentiableFunctionSetPtr_t stack =
    DifferentiableFunctionSet::create("Stack");
  stack->add (Position::create ("Position", device, ee1, MId, MId));
  stack->add (RelativeOrientation::create (
        "RelativeOrientation", device, ee1, ee2, MId,
        list_of(false)(true)(true).convert_to_container<BoolVector_t>()));
  functions.push_back (stack);
  //*/

  std::vector<ConfigurationPtr_t> cfgs (NUMBER_JACOBIAN_CALCULUS);
  for (size_t i = 0; i < NUMBER_JACOBIAN_CALCULUS; i++) cfgs[i] = cs.shoot();
  Configuration_t q1(device->currentConfiguration());
  matrix_t jacobian, fdCentral, fdForward, errorJacobian;
  for (DFs::iterator fit = functions.begin(); fit != functions.end(); ++fit) {
    DifferentiableFunction& f = **fit;
    jacobian.resize(f.outputDerivativeSize (), f.inputDerivativeSize ());
    fdForward.resize(f.outputDerivativeSize (), f.inputDerivativeSize ());
    fdCentral.resize(f.outputDerivativeSize (), f.inputDerivativeSize ());

    for (size_t i = 0; i < NUMBER_JACOBIAN_CALCULUS; i++) {
      q1 = *cfgs[i];
      jacobian.setZero ();
      f.jacobian (jacobian, q1);

      const value_type eps = std::sqrt(Eigen::NumTraits<value_type>::epsilon());

      // fdForward.setZero(); f.finiteDifferenceForward(fdForward, q1, device, eps);
      fdCentral.setZero(); f.finiteDifferenceCentral(fdCentral, q1, device, eps);

      // Forward: check the error
      // errorJacobian = jacobian - fdForward;
      // checkJacobianDiffIsZero<true> (f.name(), errorJacobian, eps);

      // Central: check the error
      errorJacobian = jacobian - fdCentral;
      checkJacobianDiffIsZero<false> (f.name(), errorJacobian, sqrt(eps));
    }
  }
}

BOOST_AUTO_TEST_CASE (SymbolicCalculus_position) {
  DevicePtr_t device = createRobot ();
  JointPtr_t ee1 = device->getJointByName ("lleg5_joint"),
             ee2 = device->getJointByName ("rleg5_joint");
  BOOST_REQUIRE (device);
  BasicConfigurationShooter cs (device);

  /// Create the constraints
  typedef DifferentiableFunctionPtr_t DFptr;
  DFptr pos = Position::create ("Position", device, ee1, MId, MId);
  Traits<PointInJoint>::Ptr_t pij  = PointInJoint::create (ee1, vector3_t(0,0,0));
  Traits<PointInJoint>::Ptr_t pij2 = PointInJoint::create (ee2, vector3_t(0,0,0));
  Traits<CalculusBaseAbstract<> >::Ptr_t relpos_sb_ptr =
    JointTranspose (ee1) * (pij2 - pij);
  DFptr relpos = RelativePosition::create ("RelPos", device, ee1, ee2, MId, MId);

  ConfigurationPtr_t q1, q2 = cs.shoot ();
  matrix_t jacobian = matrix_t (pos->outputSize (), device->numberDof ());
  for (int i = 0; i < 100; i++) {
      q1 = cs.shoot ();
      device->currentConfiguration (*q1);
      device->computeForwardKinematics ();

      pij->invalidate ();
      relpos_sb_ptr->invalidate ();

      /// Position
      LiegroupElement value = (*pos) (*q1);
      pij->computeValue ();
      BOOST_CHECK (pij->value ().isApprox (value.vector()));
      jacobian.setZero ();
      pos->jacobian (jacobian, *q1);
      pij->computeJacobian ();
      BOOST_CHECK (pij->jacobian ().isApprox (jacobian));
      // Relative position
      value = (*relpos) (*q1);
      relpos_sb_ptr->computeValue ();
      BOOST_CHECK (relpos_sb_ptr->value ().isApprox (value.vector()));
      jacobian.setZero ();
      relpos->jacobian (jacobian, *q1);
      relpos_sb_ptr->computeJacobian ();
      BOOST_CHECK (relpos_sb_ptr->jacobian ().isApprox (jacobian));
  }
}

BOOST_AUTO_TEST_CASE (SymbolicCalculus_jointframe) {
  DevicePtr_t device = createRobot ();
  JointPtr_t ee1 = device->getJointByName ("lleg5_joint"),
             ee2 = device->getJointByName ("rleg5_joint");
  BOOST_REQUIRE (device);
  BasicConfigurationShooter cs (device);

  /// Create the constraints
  typedef DifferentiableFunctionPtr_t DFptr;
  DFptr trans = Transformation::create ("Transform", device, ee1, MId);
  Traits<JointFrame>::Ptr_t jf  = JointFrame::create (ee1);
  DFptr sf = SymbolicFunction<JointFrame>::create ("SymbolicFunctionTest", device, jf);

  ConfigurationPtr_t q1, q2 = cs.shoot ();
  matrix_t jacobian1 = matrix_t (trans->outputSize (), device->numberDof ());
  matrix_t jacobian2 = matrix_t (trans->outputSize (), device->numberDof ());
  for (int i = 0; i < 100; i++) {
      q1 = cs.shoot ();
      device->currentConfiguration (*q1);
      device->computeForwardKinematics ();

      LiegroupElement value1 = (*trans) (*q1);
      LiegroupElement value2 = (*sf) (*q1);
      BOOST_CHECK (value1.vector().isApprox ( value2.vector()));
      jacobian1.setZero ();
      jacobian2.setZero ();
      trans->jacobian (jacobian1, *q1);
      sf->jacobian (jacobian2, *q1);
      BOOST_CHECK (jacobian1.isApprox ( jacobian2));
  }
}
