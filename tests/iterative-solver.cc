// Copyright (c) 2017, Joseph Mirabel
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

#define BOOST_TEST_MODULE HIERARCHICAL_ITERATIVE_SOLVER
#include <boost/test/unit_test.hpp>

#include <hpp/constraints/iterative-solver.hh>

#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/simple-device.hh>
#include <hpp/constraints/generic-transformation.hh>

using namespace hpp::constraints;

BOOST_AUTO_TEST_CASE(one_layer)
{
  DevicePtr_t device = hpp::pinocchio::humanoidSimple ("test");
  BOOST_REQUIRE (device);
  JointPtr_t ee1 = device->getJointByName ("lleg5_joint"),
             ee2 = device->getJointByName ("rleg5_joint");

  Configuration_t q = device->currentConfiguration (),
                  qrand = se3::randomConfiguration(device->model());

  HierarchicalIterativeSolver solver;
  solver.maxIterations(20);
  solver.errorThreshold(1e-3);
  solver.integration(boost::bind(hpp::pinocchio::integrate<true, se3::LieGroupTpl>, device, _1, _2, _3));

  device->currentConfiguration (qrand);
  device->computeForwardKinematics ();
  Transform3f tf1 (ee1->currentTransformation ());
  Transform3f tf2 (ee2->currentTransformation ());

  solver.addStack();
  BOOST_CHECK(solver.numberStacks() == 1);

  DifferentiableFunctionStack& stack = solver.stack(0);
  stack.add(Orientation::create ("Orientation", device, ee2, tf2));
  stack.add(Position::create    ("Position"   , device, ee2, tf2, tf1));

  solver.update();

  BOOST_CHECK(solver.solve(q));
}

