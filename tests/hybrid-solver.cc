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

#define BOOST_TEST_MODULE HYBRID_SOLVER
#include <boost/test/unit_test.hpp>

#include <hpp/constraints/hybrid-solver.hh>

#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/simple-device.hh>

#include <hpp/constraints/affine-function.hh>
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

    ExplicitSolver::ColBlockIndexes inDer () const
    {
      ExplicitSolver::ColBlockIndexes ret;
      return ret;
    }

    ExplicitSolver::RowBlockIndexes outDer () const
    {
      ExplicitSolver::RowBlockIndexes ret;
      ret.addRow (idx_ - 1, length_);
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
    R += 2 * std::pow(sin(theta/2),2) / std::pow(theta,2) * xCross * xCross;
  }
  return R;
}

class ExplicitTransformation : public DifferentiableFunction
{
  public:
    JointPtr_t joint_;
    size_type in_, inDer_;
    RelativeTransformationPtr_t rt_;

    ExplicitTransformation(JointPtr_t joint, size_type in, size_type l, size_type inDer, size_type lDer)
      : DifferentiableFunction(l, lDer, 7, 6, "ExplicitTransformation"),
        joint_ (joint), in_ (in), inDer_ (inDer)
    {
      rt_ = RelativeTransformation::create("RT", joint_->robot(),
          joint_->robot()->rootJoint(),
          joint_,
          Transform3f::Identity());
    }

    ExplicitSolver::RowBlockIndexes inArg () const
    {
      ExplicitSolver::RowBlockIndexes ret;
      ret.addRow(in_, inputSize());
      return ret;
    }

    ExplicitSolver::RowBlockIndexes outArg () const
    {
      ExplicitSolver::RowBlockIndexes ret;
      ret.addRow (0, 7);
      return ret;
    }

    ExplicitSolver::ColBlockIndexes inDer () const
    {
      ExplicitSolver::ColBlockIndexes ret;
      ret.addCol(inDer_, inputDerivativeSize());
      return ret;
    }

    ExplicitSolver::RowBlockIndexes outDer () const
    {
      ExplicitSolver::RowBlockIndexes ret;
      ret.addRow (0, 6);
      return ret;
    }

    vector_t config (vectorIn_t arg) const
    {
      vector_t q = joint_->robot()->neutralConfiguration();
      q.segment(in_, inputSize()) = arg;
      return q;
      // joint_->robot()->currentConfiguration(q);
      // joint_->robot()->computeForwardKinematics();
    }

    void impl_compute (vectorOut_t result,
                       vectorIn_t arg) const
    {
      // forwardKinematics(arg);
      vector_t value(6);
      vector_t q = config(arg);
      rt_->value(value, q);
      result.head<3>() = value.head<3>();
      result.tail<4>() = Eigen::Quaternion<value_type>(exponential(value.tail<3>())).coeffs();

      // Transform3f tf1 = joint_->robot()->rootJoint()->currentTransformation();
      // Transform3f tf2 = joint_->currentTransformation();
      // Transform3f tf = tf2.inverse() * tf1;

      // result.head<3> = tf.translation();
      // result.tail<4> = Eigen::Quaternion<value_type>(tf.rotation());
    }

    void impl_jacobian (matrixOut_t jacobian,
                        vectorIn_t arg) const
    {
      // forwardKinematics(arg);
      matrix_t J(6, rt_->inputDerivativeSize());
      vector_t q = config(arg);
      rt_->jacobian(J, q);

      inDer().rview(J).writeTo(jacobian);
    }
};

typedef boost::shared_ptr<LockedJoint> LockedJointPtr_t;
typedef boost::shared_ptr<ExplicitTransformation> ExplicitTransformationPtr_t;

BOOST_AUTO_TEST_CASE(functions)
{
  HybridSolver solver(3, 3);

  /// System:
  /// f (q1, q3) = 0
  /// q2 = g(q3)
  Eigen::Matrix<value_type, 2, 3> Jf; Jf << 1, 0, 0, 0, 0, 1;
  solver.add(AffineFunctionPtr_t(new AffineFunction (Jf)), 0);

  Eigen::Matrix<value_type,1,1> Jg; Jg (0,0) = 1;
  Eigen::RowBlockIndexes inArg; inArg.addRow (2,1);
  Eigen::ColBlockIndexes inDer; inDer.addCol (2,1);
  Eigen::RowBlockIndexes outArg; outArg.addRow (0,1);
  solver.explicitSolver().add(AffineFunctionPtr_t(new AffineFunction (Jg)),
      inArg, outArg, inDer, outArg);

  solver.explicitSolverHasChanged();
  BOOST_CHECK_EQUAL(solver.dimension(), 2);

  // We add to the system h(q3) = 0
  // This function should not be removed from the system.
  Eigen::Matrix<value_type, 1, 3> Jh; Jh << 0, 0, 1;
  solver.add(AffineFunctionPtr_t(new AffineFunction (Jh)), 0);
  BOOST_CHECK_EQUAL(solver.dimension(), 3);
}

BOOST_AUTO_TEST_CASE(hybrid_solver)
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
             ee2 = device->getJointByName ("RAnkleRoll"),
             ee3 = device->getJointByName ("LWristPitch");

  Configuration_t q = device->currentConfiguration (),
                  qrand = se3::randomConfiguration(device->model());

  HybridSolver solver(device->configSize(), device->numberDof());
  solver.maxIterations(20);
  solver.errorThreshold(1e-3);
  solver.integration(boost::bind(hpp::pinocchio::integrate<true, se3::LieGroupTpl>, device, _1, _2, _3));

  device->currentConfiguration (q);
  device->computeForwardKinematics ();
  Transform3f tf1 (ee1->currentTransformation ());
  Transform3f tf2 (ee2->currentTransformation ());
  Transform3f tf3 (ee3->currentTransformation ());

  solver.add(Orientation::create ("Orientation RAnkleRoll" , device, ee2, tf2), 0);
  solver.add(Orientation::create ("Orientation LWristPitch", device, ee3, tf3), 0);
  // solver.add(Position::create    ("Position"   , device, ee2, tf2), 0);

  BOOST_CHECK(solver.numberStacks() == 1);

  ExplicitTransformationPtr_t et;
  {
    // Find a joint such that the config parameters for the chain from the root
    // joint to it are the n first parameters (i.e. q.segment(0, n)).
    // We take the one which gives the longest block
    JointPtr_t parent = device->rootJoint(), current = device->getJointAtConfigRank(7);
    while (current->parentJoint()->index() == parent->index()) {
      parent = current;
      current = device->getJointAtConfigRank(current->rankInConfiguration() + current->configSize());
    }
    // std::cout << parent->name() << std::endl;

    et.reset (new ExplicitTransformation (parent, 7, 6,
          parent->rankInConfiguration() + parent->configSize() - 7,
          parent->rankInVelocity()      + parent->numberDof () - 6));
  }

  BOOST_CHECK(solver.explicitSolver().add (et, et->inArg(), et->outArg(), et->inDer(), et->outDer()));
  solver.explicitSolverHasChanged();

  // BOOST_CHECK_EQUAL(solver.solve<lineSearch::Backtracking  >(q), HybridSolver::SUCCESS);

  Configuration_t tmp = qrand;
  BOOST_CHECK_EQUAL(solver.solve<lineSearch::Backtracking  >(qrand), HybridSolver::SUCCESS);
  qrand = tmp;
  BOOST_CHECK_EQUAL(solver.solve<lineSearch::ErrorNormBased>(qrand), HybridSolver::SUCCESS);
  qrand = tmp;
  BOOST_CHECK_EQUAL(solver.solve<lineSearch::FixedSequence >(qrand), HybridSolver::SUCCESS);

  vector_t dq (device->numberDof());
  dq.setRandom();
  qrand = tmp;
  solver.projectOnKernel (qrand, dq, tmp);
}

