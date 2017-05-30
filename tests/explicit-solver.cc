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

#define BOOST_TEST_MODULE EXPLICIT_SOLVER
#include <boost/test/unit_test.hpp>

#include <hpp/constraints/explicit-solver.hh>

#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/simple-device.hh>
#include <hpp/constraints/generic-transformation.hh>

using namespace hpp::constraints;

class LockedJoint : public DifferentiableFunction
{
  public:
    size_type idx_, length_;
    vector_t value_;

    LockedJoint(size_type idx, size_type length, vector_t value)
      : DifferentiableFunction(0, 0, length, length, "LockedJoint"),
        idx_ (idx), length_ (length), value_ (value)
    {}

    ExplicitSolver::RowBlockIndexes inArg () const
    {
      ExplicitSolver::RowBlockIndexes ret;
      return ret;
    }

    ExplicitSolver::RowBlockIndexes outArg () const
    {
      ExplicitSolver::RowBlockIndexes ret;
      ret.addRow (idx_, length_);
      return ret;
    }

    void impl_compute (vectorOut_t result,
                       vectorIn_t ) const
    {
      result = value_;
    }

    void impl_jacobian (matrixOut_t,
                        vectorIn_t ) const
    {
      // jacobian.setIdentity();
    }
};

class TestFunction : public DifferentiableFunction
{
  public:
    size_type idxIn_, idxOut_, length_;

    TestFunction(size_type idxIn, size_type idxOut, size_type length)
      : DifferentiableFunction(length, length, length, length, "TestFunction"),
        idxIn_ (idxIn), idxOut_ (idxOut), length_ (length)
    {}

    ExplicitSolver::RowBlockIndexes inArg () const
    {
      ExplicitSolver::RowBlockIndexes ret;
      ret.addRow(idxIn_, length_);
      return ret;
    }

    ExplicitSolver::RowBlockIndexes outArg () const
    {
      ExplicitSolver::RowBlockIndexes ret;
      ret.addRow (idxOut_, length_);
      return ret;
    }

    void impl_compute (vectorOut_t result,
                       vectorIn_t arg) const
    {
      result = arg;
    }

    void impl_jacobian (matrixOut_t jacobian,
                        vectorIn_t) const
    {
      jacobian.setIdentity();
    }
};

typedef boost::shared_ptr<LockedJoint> LockedJointPtr_t;
typedef boost::shared_ptr<TestFunction> TestFunctionPtr_t;

BOOST_AUTO_TEST_CASE(locked_joints)
{
  DevicePtr_t device = hpp::pinocchio::unittest::makeDevice (hpp::pinocchio::unittest::HumanoidRomeo);
  BOOST_REQUIRE (device);
  device->rootJoint()->lowerBound (0, -1);
  device->rootJoint()->lowerBound (1, -1);
  device->rootJoint()->lowerBound (2, -1);
  device->rootJoint()->upperBound (0,  1);
  device->rootJoint()->upperBound (1,  1);
  device->rootJoint()->upperBound (2,  1);

  JointPtr_t ee1 = device->getJointByName ("LAnkleRoll"),
             ee2 = device->getJointByName ("RAnkleRoll");

  LockedJointPtr_t l1 (new LockedJoint (ee1->rankInConfiguration(), 1, vector_t::Zero(1)));
  LockedJointPtr_t l2 (new LockedJoint (ee2->rankInConfiguration(), 1, vector_t::Zero(1)));
  TestFunctionPtr_t t1 (new TestFunction (ee1->rankInConfiguration(), ee2->rankInConfiguration(), 1));

  Configuration_t q = device->currentConfiguration (),
                  qrand = se3::randomConfiguration(device->model());

  {
    ExplicitSolver solver (device->configSize(), device->numberDof());
    BOOST_CHECK( solver.add(l1, l1->inArg(), l1->outArg()));
    BOOST_CHECK(!solver.add(l1, l1->inArg(), l1->outArg()));
    BOOST_CHECK( solver.add(l2, l2->inArg(), l2->outArg()));

    BOOST_CHECK(solver.solve(qrand));
    BOOST_CHECK_EQUAL(qrand[ee1->rankInConfiguration()], 0);
    BOOST_CHECK_EQUAL(qrand[ee2->rankInConfiguration()], 0);
  }

  {
    ExplicitSolver solver (device->configSize(), device->numberDof());
    BOOST_CHECK( solver.add(l1, l1->inArg(), l1->outArg()));
    BOOST_CHECK( solver.add(t1, t1->inArg(), t1->outArg()));

    BOOST_CHECK(solver.solve(qrand));
    BOOST_CHECK_EQUAL(qrand[ee1->rankInConfiguration()], 0);
    BOOST_CHECK_EQUAL(qrand[ee2->rankInConfiguration()], 0);
  }
}


