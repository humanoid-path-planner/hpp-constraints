// Copyright (c) 2017, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#define BOOST_TEST_MODULE EXPLICIT_CONSTRAINT_SET
#include <boost/test/unit_test.hpp>

#include <hpp/constraints/explicit-constraint-set.hh>
#include <hpp/constraints/explicit.hh>
#include <hpp/constraints/explicit/relative-pose.hh>
#include <hpp/constraints/locked-joint.hh>

#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/liegroup.hh>
#include <hpp/pinocchio/liegroup-element.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/simple-device.hh>
#include <hpp/pinocchio/urdf/util.hh>

#include <hpp/constraints/affine-function.hh>
#include <hpp/constraints/explicit-constraint-set.hh>
#include <hpp/constraints/generic-transformation.hh>
#include <hpp/constraints/symbolic-calculus.hh>

#include <../tests/util.hh>

using Eigen::RowBlockIndices;
using Eigen::ColBlockIndices;
using Eigen::BlockIndex;

using hpp::pinocchio::unittest::HumanoidSimple;
using hpp::pinocchio::unittest::makeDevice;
using hpp::constraints::AffineFunction;
using hpp::constraints::size_type;
using hpp::constraints::vector_t;
using hpp::constraints::matrix_t;
using hpp::constraints::ExplicitConstraintSet;
using hpp::constraints::Explicit;
using hpp::constraints::ExplicitPtr_t;
using hpp::constraints::matrix3_t;
using hpp::constraints::vector3_t;
using hpp::constraints::DifferentiableFunction;
using hpp::constraints::JointPtr_t;
using hpp::constraints::RelativeTransformation;
using hpp::constraints::RelativeTransformationPtr_t;
using hpp::constraints::LiegroupSpace;
using hpp::constraints::LiegroupSpacePtr_t;
using hpp::constraints::Transform3f;
using hpp::constraints::vectorIn_t;
using hpp::constraints::LiegroupElement;
using hpp::constraints::LiegroupElementRef;
using hpp::constraints::matrixOut_t;
using hpp::constraints::AffineFunctionPtr_t;
using hpp::constraints::segment_t;
using hpp::constraints::segments_t;
using hpp::constraints::DevicePtr_t;
using hpp::constraints::Device;
using hpp::constraints::Configuration_t;
using hpp::constraints::value_type;
using hpp::constraints::Equality;
using hpp::constraints::ComparisonTypes_t;
using hpp::constraints::LockedJoint;
using hpp::constraints::LockedJointPtr_t;
using hpp::constraints::EqualToZero;
using hpp::constraints::Equality;

namespace Eigen {
  namespace internal {
    bool operator== (const empty_struct&, const empty_struct&) { return true; }
  } // namespace internal

  template <bool _allRows, bool _allCols>
  bool operator== (const MatrixBlocks<_allRows,_allCols>& a,
                   const MatrixBlocks<_allRows,_allCols>& b)
  {
    return ( _allRows || a.nbRows() == b.nbRows())
      &&   ( _allCols || a.nbCols() == b.nbCols())
      &&   ( _allRows || a.rows()   == b.rows())
      &&   ( _allCols || a.cols()   == b.cols());
  }
} // namespace Eigen

// Explicit expression of variables
//   - [idxOut:idxOut+length] with respect to
//   - [idxIn:idxIn+length]
class TestFunction : public AffineFunction
{
  public:
    size_type idxIn_, idxOut_, length_;

    TestFunction(size_type idxIn, size_type idxOut, size_type length)
      : AffineFunction(matrix_t::Identity(length,length), "TestFunction"),
        idxIn_ (idxIn), idxOut_ (idxOut), length_ (length)
    {}

    // return configuration input interval
    ExplicitConstraintSet::RowBlockIndices inArg () const
    {
      ExplicitConstraintSet::RowBlockIndices ret;
      ret.addRow(idxIn_, length_);
      return ret;
    }

    // return configuration output interval
    ExplicitConstraintSet::RowBlockIndices outArg () const
    {
      ExplicitConstraintSet::RowBlockIndices ret;
      ret.addRow (idxOut_, length_);
      return ret;
    }

    // return velocity input interval
    ExplicitConstraintSet::ColBlockIndices inDer () const
    {
      ExplicitConstraintSet::ColBlockIndices ret;
      ret.addCol(idxIn_ - 1, length_); // TODO this assumes there is only the freeflyer
      return ret;
    }

    // return velocity output interval
    ExplicitConstraintSet::RowBlockIndices outDer () const
    {
      ExplicitConstraintSet::RowBlockIndices ret;
      ret.addRow (idxOut_ - 1, length_); // TODO this assumes there is only the freeflyer
      return ret;
    }
}; // class TestFunction

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

    ExplicitTransformation(JointPtr_t joint, size_type in, size_type l,
                           size_type inDer, size_type lDer)
      : DifferentiableFunction(l, lDer,
                               LiegroupSpacePtr_t (LiegroupSpace::R3xSO3 ()),
                               "ExplicitTransformation"),
        joint_ (joint), in_ (in), inDer_ (inDer)
    {
      rt_ = RelativeTransformation::create("RT", joint_->robot(),
          joint_->robot()->rootJoint(),
          joint_,
          Transform3f::Identity());
    }

    ExplicitConstraintSet::RowBlockIndices inArg () const
    {
      ExplicitConstraintSet::RowBlockIndices ret;
      ret.addRow(in_, inputSize());
      return ret;
    }

    ExplicitConstraintSet::RowBlockIndices outArg () const
    {
      ExplicitConstraintSet::RowBlockIndices ret;
      ret.addRow (0, 7);
      return ret;
    }

    ExplicitConstraintSet::ColBlockIndices inDer () const
    {
      ExplicitConstraintSet::ColBlockIndices ret;
      ret.addCol(inDer_, inputDerivativeSize());
      return ret;
    }

    ExplicitConstraintSet::RowBlockIndices outDer () const
    {
      ExplicitConstraintSet::RowBlockIndices ret;
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

    void impl_compute (LiegroupElementRef result,
                       vectorIn_t arg) const
    {
      // forwardKinematics(arg);
      LiegroupElement transform (LiegroupSpace::Rn (6));
      vector_t q = config(arg);
      rt_->value (transform, q);
      result. vector ().head<3>() = transform. vector ().head<3>();
      result. vector ().tail<4>() =
        Eigen::Quaternion<value_type>
        (exponential(transform. vector ().tail<3>())).coeffs();
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

typedef hpp::shared_ptr<TestFunction> TestFunctionPtr_t;
typedef hpp::shared_ptr<ExplicitTransformation> ExplicitTransformationPtr_t;

template <int N>
void order_test (const AffineFunctionPtr_t f[N], const segment_t s[N+1],
    const std::vector<int> forder,
    const segments_t& inArgs,
    const segments_t& outArgs)
{
  ExplicitConstraintSet expression (LiegroupSpace::Rn (4));
  for (int i = 0; i < N; ++i) {
    int fo = forder[i],
    si = forder[i], so = forder[i] + 1;
    ExplicitPtr_t constraint (Explicit::create
                              (LiegroupSpace::Rn (4), f[fo],
                               segments_t (1, s[si]),
                               segments_t (1, s[so]),
                               segments_t (1, s[si]),
                               segments_t (1, s[so])));
    BOOST_CHECK( expression.add(constraint) >= 0);
  }
  BOOST_CHECK_EQUAL( expression.inArgs().rows(), inArgs);
  BOOST_CHECK_EQUAL( expression.outArgs().rows(), outArgs);
}

BOOST_AUTO_TEST_CASE(order)
{
  Eigen::Matrix<value_type,1,1> M; M(0,0) = 1;

  // dof     :  0 -> 1 -> 2 -> 3
  // function:    f0   f1   f2
  AffineFunctionPtr_t f[] = {
    AffineFunction::create (M),
    AffineFunction::create (M),
    AffineFunction::create (M)
  };
  segment_t s[] = { segment_t (0, 1), segment_t (1, 1), segment_t (2, 1), segment_t (3, 1) };
  segments_t inArgs { s[0] },
             outArgs { s[1], s[2], s[3] };
  BlockIndex::shrink (outArgs);

  std::vector<int> order(3);

  order = { 0, 1, 2 };
  order_test<3> (f, s, order, inArgs, outArgs);
  order = { 0, 2, 1 };
  order_test<3> (f, s, order, inArgs, outArgs);
  order = { 1, 0, 2 };
  order_test<3> (f, s, order, inArgs, outArgs);
  order = { 1, 2, 0 };
  order_test<3> (f, s, order, inArgs, outArgs);
  order = { 2, 0, 1 };
  order_test<3> (f, s, order, inArgs, outArgs);
  order = { 2, 1, 0 };
  order_test<3> (f, s, order, inArgs, outArgs);
}

BOOST_AUTO_TEST_CASE(jacobian1)
{
  matrix_t J[] = {
      (matrix_t(1,1) << 1).finished()
    , (matrix_t(1,1) << 2).finished()
    , (matrix_t(1,1) << 3).finished()
  };

  // dof     :  0 -> 1 -> 2 -> 3
  // function:    f0   f1   f2
  AffineFunctionPtr_t f[] = {
    AffineFunction::create(J[0]),
    AffineFunction::create(J[1]),
    AffineFunction::create(J[2])
  };
  segments_t s[] = { segments_t (1, segment_t (0, 1)),
                     segments_t (1, segment_t (1, 1)),
                     segments_t (1, segment_t (2, 1)),
                     segments_t (1, segment_t (3, 1)) };

  ExplicitConstraintSet expression (LiegroupSpace::Rn (4));
  for (int i = 0; i < 3; ++i) {
    ExplicitPtr_t constraint (Explicit::create
                              (LiegroupSpace::Rn (4), f[i], s[i], s[i+1], s[i],
                               s[i+1]));
    expression.add (constraint);
  }

  vector_t x(4); x << 1,2,3,4;
  vector_t xres = x;
  BOOST_CHECK (expression.solve(xres));

  // Check the solution
  BOOST_CHECK_EQUAL (xres[0], x[0]);
  for (int i = 0; i < 3; ++i)
    BOOST_CHECK_EQUAL (xres.segment<1>(i+1), (*f[i])(xres.segment<1>(i)).vector());

  // Check the jacobian
  // It should be ( J[0], J[1] * J[0], J[2] * J[1] * J[0])
  matrix_t expjac (matrix_t::Zero(expression.nv (), expression.nv ()));
  expjac.col(0) << 1, J[0], J[1] * J[0], J[2] * J[1] * J[0];
  matrix_t jacobian (expression.nv (), expression.nv ());
  expression.jacobian (jacobian, xres);
  BOOST_CHECK_EQUAL (jacobian, expjac);
}

BOOST_AUTO_TEST_CASE(jacobian2)
{
  matrix_t J[] = {
      (matrix_t(1,2) << 3.2, -0.3).finished()
    , (matrix_t(1,1) << 4.1).finished()
    , (matrix_t(1,2) << -0.3, 1.2).finished()
  };

  /* dof     :  1,2 -> 0 \
   * function:      f0    --> 4
   *            1   -> 3 / f2
   *                f1
   */
  AffineFunctionPtr_t f[] = {
    AffineFunction::create (J[0])
    , AffineFunction::create(J[1])
    , AffineFunction::create(J[2])
  };
  std::vector<segments_t> s(6);
  s[0] = { segment_t (1, 2) };
  s[1] = { segment_t (0, 1) };
  s[2] = { segment_t (3, 1) };
  s[3] = { segment_t (4, 1) };
  s[4] = { segment_t (0, 1), segment_t (3, 1) };
  s[5] = { segment_t (1, 1) };

  ExplicitConstraintSet expression (LiegroupSpace::Rn (5));
  ExplicitPtr_t constraint;
  constraint = Explicit::create
    (LiegroupSpace::Rn (5), f[0], s[0], s[1], s[0], s[1]);
  expression.add (constraint);
  constraint = Explicit::create
    (LiegroupSpace::Rn (5), f[2], s[4], s[3], s[4], s[3]);
  expression.add (constraint);
  constraint = Explicit::create
    (LiegroupSpace::Rn (5), f[1], s[5], s[2], s[5], s[2]);
  expression.add (constraint);

  Eigen::MatrixXi inOutDependencies (3, 5);
  inOutDependencies << 0, 1, 1, 0, 0,
                       0, 2, 1, 0, 0,
                       0, 1, 0, 0, 0;
  BOOST_CHECK_EQUAL (expression.inOutDependencies(), inOutDependencies);
  inOutDependencies.resize (3, 2);
  inOutDependencies << 1, 1,
                       1, 0,
                       2, 1;
  BOOST_CHECK_EQUAL (expression.inOutDofDependencies(), inOutDependencies);

  segments_t inArgs = s[0],
             outArgs = { s[1][0], s[2][0], s[3][0] };
  BlockIndex::shrink (outArgs);

  BOOST_CHECK_EQUAL( expression.inArgs().rows(), inArgs);
  BOOST_CHECK_EQUAL( expression.outArgs().rows(), outArgs);

  vector_t x(5); x << 1,2,3,4,5;
  vector_t xres = x;
  BOOST_CHECK (expression.solve(xres));

  // Check the solution
  BOOST_CHECK_EQUAL (xres.segment<2>(1), x.segment<2>(1));
  BOOST_CHECK_EQUAL (xres.segment<1>(0), (*f[0])(xres.segment<2>(1)).vector());
  BOOST_CHECK_EQUAL (xres.segment<1>(3), (*f[1])(xres.segment<1>(1)).vector());
  BOOST_CHECK_EQUAL (xres.segment<1>(4), (*f[2])(RowBlockIndices(s[4]).rview(xres).eval()).vector());

  // Check the jacobian
  // It should be ( J[0], J[1] * J[0], J[2] * J[1] * J[0])
  matrix_t expjac (matrix_t::Zero(expression.nv (), expression.nv ()));
  expjac.block<5, 2>(0,1) <<
    J[0](0,0), J[0](0,1),
    1, 0,
    0, 1,
    J[1](0,0), 0,
    J[2](0,0) * J[0](0,0) + J[2](0,1) * J[1](0,0), J[2](0,0) * J[0](0,1);
  matrix_t jacobian (expression.nv (), expression.nv ());
  expression.jacobian (jacobian, xres);
  BOOST_CHECK_EQUAL (jacobian, expjac);

  matrix_t smallJ = expression.jacobianNotOutToOut (jacobian);
  BOOST_CHECK_EQUAL (expression.outArgs().rview(xres).eval(),
      smallJ * expression.inArgs().rview(xres).eval());
}

BOOST_AUTO_TEST_CASE(locked_joints)
{
  DevicePtr_t device (makeDevice (HumanoidSimple));
  device->controlComputation((hpp::pinocchio::Computation_t) (hpp::pinocchio::JOINT_POSITION | hpp::pinocchio::JACOBIAN));

  BOOST_REQUIRE (device);
  device->rootJoint()->lowerBound (0, -1);
  device->rootJoint()->lowerBound (1, -1);
  device->rootJoint()->lowerBound (2, -1);
  device->rootJoint()->upperBound (0,  1);
  device->rootJoint()->upperBound (1,  1);
  device->rootJoint()->upperBound (2,  1);

  JointPtr_t ee1 = device->getJointByName ("lleg5_joint"),
             ee2 = device->getJointByName ("rleg5_joint"),
             ee3 = device->getJointByName ("rleg4_joint");

  LockedJointPtr_t l1 (LockedJoint::create
                       (ee1, ee1->configurationSpace ()->neutral ()));
  LockedJointPtr_t l2 (LockedJoint::create
                       (ee2, ee2->configurationSpace ()->neutral ()));
  LockedJointPtr_t l3 (LockedJoint::create
                       (ee3, ee3->configurationSpace ()->neutral ()));
  TestFunctionPtr_t t1 (new TestFunction (ee1->rankInConfiguration(),
                                          ee2->rankInConfiguration(), 1));
  TestFunctionPtr_t t2 (new TestFunction (ee2->rankInConfiguration(),
                                          ee1->rankInConfiguration(), 1));

  RowBlockIndices expectedRow;
  ColBlockIndices expectedCol;

  Configuration_t q = device->currentConfiguration (),
                  qrand = ::pinocchio::randomConfiguration(device->model());

  {
    ExplicitConstraintSet expression (device->configSpace ());
    ExplicitPtr_t constraint;
    BOOST_CHECK (expression.add (l1) >= 0);
    BOOST_CHECK (expression.add (l1) <  0);
    BOOST_CHECK (expression.add (l2) >= 0);

    expectedRow = RowBlockIndices();
    expectedRow.addRow (ee1->rankInConfiguration(), 1);
    expectedRow.addRow (ee2->rankInConfiguration(), 1);
    expectedRow.updateRows<true,true,true>();
    BOOST_CHECK_EQUAL (expression.outArgs(), expectedRow);

    expectedRow = RowBlockIndices(BlockIndex::difference (
          BlockIndex::segment_t(0, expression.nq ()),
          expectedRow.rows()));
    BOOST_CHECK_EQUAL (expression.notOutArgs (), expectedRow);

    expectedRow = RowBlockIndices();
    expectedRow.addRow (ee1->rankInVelocity(), 1);
    expectedRow.addRow (ee2->rankInVelocity(), 1);
    expectedRow.updateRows<true,true,true>();
    BOOST_CHECK_EQUAL (expression.outDers(), expectedRow);

    expectedCol = ColBlockIndices(BlockIndex::difference (
          BlockIndex::segment_t(0, expression.nv ()),
          expectedRow.rows()));
    BOOST_CHECK_EQUAL (expression.notOutDers(), expectedCol);

    expectedRow = RowBlockIndices();
    BOOST_CHECK_EQUAL (expression.inArgs(), expectedRow);
    expectedCol = ColBlockIndices();
    BOOST_CHECK_EQUAL (expression.inDers(), expectedCol);

    BOOST_CHECK(expression.solve(qrand));
    BOOST_CHECK_EQUAL(qrand[ee1->rankInConfiguration()], 0);
    BOOST_CHECK_EQUAL(qrand[ee2->rankInConfiguration()], 0);

    expression.rightHandSide(l1, vector_t::Ones(1));
    expression.rightHandSide(l2, vector_t::Constant(1,-0.2));
    BOOST_CHECK(expression.solve(qrand));
    BOOST_CHECK_EQUAL(qrand[ee1->rankInConfiguration()], 1);
    BOOST_CHECK_EQUAL(qrand[ee2->rankInConfiguration()], -0.2);

    matrix_t jacobian (device->numberDof(), device->numberDof());
    expression.jacobian(jacobian, q);
    BOOST_CHECK(expression.jacobianNotOutToOut (jacobian).eval ().isZero ());
  }

  {
    ExplicitConstraintSet expression (device->configSpace());
    BOOST_CHECK( expression.add (l1) >= 0);
    ExplicitPtr_t constraint;
    constraint = Explicit::create
      (device->configSpace (), t1, t1->inArg().indices (),
       t1->outArg().indices (), t1->inDer().indices (),
       t1->outDer().indices ());
    BOOST_CHECK (expression.add (constraint) >= 0);

    BOOST_CHECK(expression.solve(qrand));
    vector_t error(expression.outDers().nbIndices());
    BOOST_CHECK(expression.isSatisfied(qrand, error));
    BOOST_CHECK_EQUAL(qrand[ee1->rankInConfiguration()], 0);
    BOOST_CHECK_EQUAL(qrand[ee2->rankInConfiguration()], 0);

    matrix_t jacobian (device->numberDof(), device->numberDof());
    expression.jacobian(jacobian, q);
    BOOST_CHECK(expression.jacobianNotOutToOut (jacobian).eval ().isZero ());
  }

  {
    ExplicitConstraintSet expression (device->configSpace ());
    ExplicitPtr_t constraint;
    constraint = Explicit::create
      (device->configSpace (), t1, t1->inArg().indices (),
       t1->outArg().indices (), t1->inDer().indices (),
       t1->outDer().indices ());
    BOOST_CHECK (expression.add (constraint) >= 0);

    matrix_t jacobian (device->numberDof(), device->numberDof());
    expression.jacobian(jacobian, q);
    BOOST_CHECK_EQUAL(jacobian(ee2->rankInVelocity(), ee1->rankInVelocity()), 1);
    BOOST_CHECK_EQUAL
      (expression.jacobianNotOutToOut (jacobian).eval ().norm (), 1);
  }

  {
    ExplicitConstraintSet expression (device->configSpace());
    ExplicitPtr_t constraint;
    constraint = Explicit::create
      (device->configSpace (), t1, t1->inArg().indices (),
       t1->outArg().indices (), t1->inDer().indices (),
       t1->outDer().indices ());
    BOOST_CHECK (expression.add (constraint) >= 0);
    constraint = Explicit::create
      (device->configSpace (), t2, t2->inArg().indices (),
       t2->outArg().indices (), t2->inDer().indices (),
       t2->outDer().indices ());
    BOOST_CHECK (expression.add (constraint) <  0);
  }

  {
    ExplicitConstraintSet expression (device->configSpace ());
    ExplicitPtr_t constraint;
    constraint = Explicit::create
      (device->configSpace (), t1, t1->inArg().indices (),
       t1->outArg().indices (), t1->inDer().indices (),
       t1->outDer().indices ());
    BOOST_CHECK (expression.add (constraint) >= 0);
    BOOST_CHECK (expression.add (l2) <  0);
    BOOST_CHECK (expression.add (l3) >= 0);

    matrix_t jacobian (device->numberDof(), device->numberDof());
    expression.jacobian(jacobian, q);
    BOOST_CHECK_EQUAL(jacobian(ee2->rankInVelocity(), ee1->rankInVelocity()), 1);
    BOOST_CHECK_EQUAL
      (expression.jacobianNotOutToOut (jacobian).eval ().norm (), 1);
  }

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

    ExplicitTransformationPtr_t et (new ExplicitTransformation (parent, 7, 6,
          parent->rankInConfiguration() + parent->configSize() - 7,
          parent->rankInVelocity()      + parent->numberDof () - 6));

    ExplicitConstraintSet expression (device->configSpace ());
    ExplicitPtr_t constraint;
    constraint = Explicit::create
      (device->configSpace (), et, et->inArg().indices (),
       et->outArg().indices (), et->inDer().indices (),
       et->outDer().indices ());
    BOOST_CHECK (expression.add (constraint) >= 0);
    BOOST_CHECK (expression.add (l2) >= 0);

    matrix_t jacobian (device->numberDof(), device->numberDof());
    expression.jacobian(jacobian, qrand);
  }
}

BOOST_AUTO_TEST_CASE(RelativePose)
{
  const std::string urdf
    ("<robot name=\"two-freeflyers\">\n"
     "  <link name=\"base_link\">\n"
     "  </link>\n"
     "</robot>");
  // Make robot with two free flyers.
  DevicePtr_t device (Device::create("two-freeflyers"));
  hpp::pinocchio::urdf::loadModelFromString(device, 0, "1", "freeflyer",
                                            urdf, "");
  hpp::pinocchio::urdf::loadModelFromString(device, 0, "2", "freeflyer",
                                            urdf, "");
  assert(device->configSize() == 14);
  assert(device->numberDof() == 12);
  // Set joint bounds
  hpp::pinocchio::JointPtr_t joint1 (device->jointAt(0));
  hpp::pinocchio::JointPtr_t joint2 (device->jointAt(1));
  vector_t low (device->configSize()); low.fill(-1);
  vector_t  up (device->configSize());  up.fill( 1);

  hpp::constraints::Transform3f frame1(pinocchio::SE3::Random());
  hpp::constraints::Transform3f frame2(pinocchio::SE3::Random());
  // explicit relative pose
  hpp::constraints::ExplicitPtr_t constraint
    (hpp::constraints::explicit_::RelativePose::create
     ("explicit-relative-pose", device, joint1, joint2, frame1, frame2,
      6 * Equality, std::vector<bool>(6,true)));
  assert(constraint->inputConf().size() == 1);
  assert(constraint->inputConf()[0].first == 0);
  assert(constraint->inputConf()[0].second == 7);

  for (size_type k=0; k<100; ++k) {
    // Pick a random configuration
    vector_t q_rand(pinocchio::randomConfiguration(device->model(), low, up));
    for (size_type i=0; i<64; ++i) {
      vector_t q (q_rand);
      size_type m = 1;
      hpp::constraints::ComparisonTypes_t comp(6 * EqualToZero);
      for (size_type j=0; j<6; ++j) {
        // m = 2^(j+1)
        if ((m & i) == 0) {
          comp [j] = Equality;
        } else {
          comp [j] = EqualToZero;
        }
        m *= 2;
      }
      constraint->comparisonType (comp);
      hpp::constraints::ExplicitConstraintSet ecs(device->configSpace());
      size_type res(ecs.add(constraint));
      BOOST_CHECK (res != -1);
      // Initialize right hand side with initial value of q -> rhs0_impl
      LiegroupElement rhs0_impl (constraint->function().outputSpace());
      rhs0_impl.vector() = ecs.rightHandSideFromInput(q);

      // Solve constraint
      ecs.solve(q);
      // Get value of function h for new value of q -> rhs2_impl
      LiegroupElement rhs2_impl (constraint->function().outputSpace());
      constraint->function().value(rhs2_impl, q);
      vector_t logRhs0_impl = log(rhs0_impl);
      vector_t logRhs2_impl = log(rhs2_impl);
      assert(logRhs0_impl.size()==6);
      assert(logRhs2_impl.size()==6);
      // For each coordinate check that
      // rhs2_impl [j] == 0 if comp[i] is EqualToZero
      // rhs2_impl [j] == rhs0_impl [j] if comp[i] is Equality
      for (size_type j=0; j<6; ++j){
        BOOST_CHECK(((comp[j] == EqualToZero) &&
                     (fabs(logRhs2_impl[j]) < 1e-10))
                    || ((comp[j] == Equality) &&
                        (fabs(logRhs2_impl[j] - logRhs0_impl[j])
                         < 1e-10)));
      }
    }
  }
}
