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

#define BOOST_TEST_MODULE EXPLICIT_CONSTRAINT_SET
#include <boost/test/unit_test.hpp>
#include <boost/assign/list_of.hpp>

#include <hpp/constraints/explicit-constraint-set.hh>

#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/liegroup.hh>
#include <hpp/pinocchio/liegroup-element.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/simple-device.hh>

#include <hpp/constraints/affine-function.hh>
#include <hpp/constraints/generic-transformation.hh>
#include <hpp/constraints/symbolic-calculus.hh>

#include <../tests/util.hh>

using boost::assign::list_of;

using Eigen::RowBlockIndices;
using Eigen::ColBlockIndices;
using Eigen::BlockIndex;

using hpp::pinocchio::unittest::HumanoidRomeo;
using hpp::pinocchio::unittest::makeDevice;
using hpp::constraints::AffineFunction;
using hpp::constraints::size_type;
using hpp::constraints::vector_t;
using hpp::constraints::matrix_t;
using hpp::constraints::ExplicitConstraintSet;
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

class LockedJoint : public AffineFunction
{
  public:
    size_type idx_, length_;
    vector_t value_;

    LockedJoint(size_type idx, size_type length, vector_t value)
      : AffineFunction(matrix_t(length,0), value, "LockedJoint"),
        idx_ (idx), length_ (length), value_ (value)
    {}

    ExplicitConstraintSet::RowBlockIndices inArg () const
    {
      ExplicitConstraintSet::RowBlockIndices ret;
      return ret;
    }

    ExplicitConstraintSet::RowBlockIndices outArg () const
    {
      ExplicitConstraintSet::RowBlockIndices ret;
      ret.addRow (idx_, length_);
      return ret;
    }

    ExplicitConstraintSet::ColBlockIndices inDer () const
    {
      ExplicitConstraintSet::ColBlockIndices ret;
      return ret;
    }

    ExplicitConstraintSet::RowBlockIndices outDer () const
    {
      ExplicitConstraintSet::RowBlockIndices ret;
      ret.addRow (idx_ - 1, length_);
      return ret;
    }
}; // class LockedJoint

class TestFunction : public AffineFunction
{
  public:
    size_type idxIn_, idxOut_, length_;

    TestFunction(size_type idxIn, size_type idxOut, size_type length)
      : AffineFunction(matrix_t::Identity(length,length), "TestFunction"),
        idxIn_ (idxIn), idxOut_ (idxOut), length_ (length)
    {}

    ExplicitConstraintSet::RowBlockIndices inArg () const
    {
      ExplicitConstraintSet::RowBlockIndices ret;
      ret.addRow(idxIn_, length_);
      return ret;
    }

    ExplicitConstraintSet::RowBlockIndices outArg () const
    {
      ExplicitConstraintSet::RowBlockIndices ret;
      ret.addRow (idxOut_, length_);
      return ret;
    }

    ExplicitConstraintSet::ColBlockIndices inDer () const
    {
      ExplicitConstraintSet::ColBlockIndices ret;
      ret.addCol(idxIn_ - 1, length_); // TODO this assumes there is only the freeflyer
      return ret;
    }

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

    void impl_compute (LiegroupElement& result,
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

typedef boost::shared_ptr<LockedJoint> LockedJointPtr_t;
typedef boost::shared_ptr<TestFunction> TestFunctionPtr_t;
typedef boost::shared_ptr<ExplicitTransformation> ExplicitTransformationPtr_t;

template <int N>
void order_test (const AffineFunctionPtr_t f[N], const segment_t s[N+1],
    const std::vector<int> forder,
    const segments_t& inArgs,
    const segments_t& outArgs)
{
  ExplicitConstraintSet expression (4, 4);
  for (int i = 0; i < N; ++i) {
    int fo = forder[i],
    si = forder[i], so = forder[i] + 1;
    BOOST_CHECK( expression.add(f[fo], s[si], s[so], s[si], s[so]) >= 0);
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
    AffineFunctionPtr_t (new AffineFunction (M)),
    AffineFunctionPtr_t (new AffineFunction (M)),
    AffineFunctionPtr_t (new AffineFunction (M))
  };
  segment_t s[] = { segment_t (0, 1), segment_t (1, 1), segment_t (2, 1), segment_t (3, 1) };
  segments_t inArgs = list_of(s[0]),
             outArgs = list_of(s[1])(s[2])(s[3]);
  BlockIndex::shrink (outArgs);

  std::vector<int> order(3);

  order = list_of(0)(1)(2);
  order_test<3> (f, s, order, inArgs, outArgs);
  order = list_of(0)(2)(1);
  order_test<3> (f, s, order, inArgs, outArgs);
  order = list_of(1)(0)(2);
  order_test<3> (f, s, order, inArgs, outArgs);
  order = list_of(1)(2)(0);
  order_test<3> (f, s, order, inArgs, outArgs);
  order = list_of(2)(0)(1);
  order_test<3> (f, s, order, inArgs, outArgs);
  order = list_of(2)(1)(0);
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
    AffineFunctionPtr_t (new AffineFunction (J[0])),
    AffineFunctionPtr_t (new AffineFunction (J[1])),
    AffineFunctionPtr_t (new AffineFunction (J[2]))
  };
  segment_t s[] = { segment_t (0, 1), segment_t (1, 1), segment_t (2, 1), segment_t (3, 1) };

  ExplicitConstraintSet expression (4, 4);
  for (int i = 0; i < 3; ++i)
    expression.add(f[i], s[i], s[i+1], s[i], s[i+1]);

  vector_t x(4); x << 1,2,3,4;
  vector_t xres = x;
  BOOST_CHECK (expression.solve(xres));

  // Check the solution
  BOOST_CHECK_EQUAL (xres[0], x[0]);
  for (int i = 0; i < 3; ++i)
    BOOST_CHECK_EQUAL (xres.segment<1>(i+1), (*f[i])(xres.segment<1>(i)).vector());

  // Check the jacobian
  // It should be ( J[0], J[1] * J[0], J[2] * J[1] * J[0])
  matrix_t expjac (matrix_t::Zero(expression.derSize(), expression.derSize()));
  expjac.col(0) << 1, J[0], J[1] * J[0], J[2] * J[1] * J[0];
  matrix_t jacobian (expression.derSize(), expression.derSize());
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
      AffineFunctionPtr_t (new AffineFunction (J[0]))
    , AffineFunctionPtr_t (new AffineFunction (J[1]))
    , AffineFunctionPtr_t (new AffineFunction (J[2]))
  };
  std::vector<segments_t> s(6);
  s[0] = (list_of(segment_t (1, 2)));
  s[1] = (list_of(segment_t (0, 1)));
  s[2] = (list_of(segment_t (3, 1)));
  s[3] = (list_of(segment_t (4, 1)));
  s[4] = (list_of(segment_t (0, 1))(segment_t (3, 1)));
  s[5] = (list_of(segment_t (1, 1)));

  ExplicitConstraintSet expression (5, 5);
  expression.add(f[0], s[0], s[1], s[0], s[1]);
  expression.add(f[2], s[4], s[3], s[4], s[3]);
  expression.add(f[1], s[5], s[2], s[5], s[2]);

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
             outArgs = list_of(s[1][0])(s[2][0])(s[3][0]);
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
  matrix_t expjac (matrix_t::Zero(expression.derSize(), expression.derSize()));
  expjac.block<5, 2>(0,1) <<
    J[0](0,0), J[0](0,1),
    1, 0,
    0, 1,
    J[1](0,0), 0,
    J[2](0,0) * J[0](0,0) + J[2](0,1) * J[1](0,0), J[2](0,0) * J[0](0,1);
  matrix_t jacobian (expression.derSize(), expression.derSize());
  expression.jacobian (jacobian, xres);
  BOOST_CHECK_EQUAL (jacobian, expjac);

  matrix_t smallJ = expression.viewJacobian(jacobian);
  BOOST_CHECK_EQUAL (expression.outArgs().rview(xres).eval(),
      smallJ * expression.inArgs().rview(xres).eval());
}

BOOST_AUTO_TEST_CASE(jacobian3)
{
  matrix_t J[] = {
      (matrix_t(1,2) << 3.2, -0.3).finished()
    , (matrix_t(1,1) << 4.1).finished()
    , (matrix_t(1,2) << -0.3, 1.2).finished()
  };
  matrix_t Jg0   (1,1); Jg0    << 2;
  matrix_t Jginv0(1,1); Jginv0 << 0.5;

  /* dof     :  1,2 -> 0 \
   * function:      f0    --> 4
   *            1   -> 3 / f2
   *                f1
   */
  AffineFunctionPtr_t f[] = {
      AffineFunctionPtr_t (new AffineFunction (J[0]))
    , AffineFunctionPtr_t (new AffineFunction (J[1]))
    , AffineFunctionPtr_t (new AffineFunction (J[2]))
  };
  AffineFunctionPtr_t g0    (new AffineFunction (Jg0   )),
                      ginv0 (new AffineFunction (Jginv0));
  std::vector<segments_t> s(6);
  s[0] = (list_of(segment_t (1, 2)));
  s[1] = (list_of(segment_t (0, 1)));
  s[2] = (list_of(segment_t (3, 1)));
  s[3] = (list_of(segment_t (4, 1)));
  s[4] = (list_of(segment_t (0, 1))(segment_t (3, 1)));
  s[5] = (list_of(segment_t (1, 1)));

  ExplicitConstraintSet expression (5, 5);
  expression.add(f[0], s[0], s[1], s[0], s[1]);
  expression.add(f[2], s[4], s[3], s[4], s[3]);
  expression.add(f[1], s[5], s[2], s[5], s[2]);
  expression.setG(f[0], g0, ginv0);

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
             outArgs = list_of(s[1][0])(s[2][0])(s[3][0]);
  BlockIndex::shrink (outArgs);

  BOOST_CHECK_EQUAL( expression.inArgs().rows(), inArgs);
  BOOST_CHECK_EQUAL( expression.outArgs().rows(), outArgs);

  vector_t x(5); x << 1,2,3,4,5;
  vector_t xres = x;
  BOOST_CHECK (expression.solve(xres));

  // Check the solution
  BOOST_CHECK_EQUAL (xres.segment<2>(1), x.segment<2>(1));
  BOOST_CHECK_EQUAL (xres.segment<1>(0), Jginv0 * (*f[0])(xres.segment<2>(1)).vector());
  BOOST_CHECK_EQUAL (xres.segment<1>(3), (*f[1])(xres.segment<1>(1)).vector());
  BOOST_CHECK_EQUAL (xres.segment<1>(4), (*f[2])(RowBlockIndices(s[4]).rview(xres).eval()).vector());

  // Check the jacobian
  // It should be ( J[0], J[1] * J[0], J[2] * J[1] * J[0])
  matrix_t expjac (matrix_t::Zero(expression.derSize(), expression.derSize()));
  expjac.block<5, 2>(0,1) <<
    Jginv0(0) * J[0](0,0), Jginv0(0) * J[0](0,1),
    1, 0,
    0, 1,
    J[1](0,0), 0,
    J[2](0,0) * Jginv0(0) * J[0](0,0) + J[2](0,1) * J[1](0,0), J[2](0,0) * Jginv0(0) * J[0](0,1);
  matrix_t jacobian (expression.derSize(), expression.derSize());
  expression.jacobian (jacobian, xres);
  BOOST_CHECK_EQUAL (jacobian, expjac);

  matrix_t smallJ = expression.viewJacobian(jacobian);
  BOOST_CHECK_EQUAL (expression.outArgs().rview(xres).eval(),
      smallJ * expression.inArgs().rview(xres).eval());
}

BOOST_AUTO_TEST_CASE(locked_joints)
{
  DevicePtr_t device (makeDevice (HumanoidRomeo));
  device->controlComputation((Device::Computation_t) (Device::JOINT_POSITION | Device::JACOBIAN));

  BOOST_REQUIRE (device);
  device->rootJoint()->lowerBound (0, -1);
  device->rootJoint()->lowerBound (1, -1);
  device->rootJoint()->lowerBound (2, -1);
  device->rootJoint()->upperBound (0,  1);
  device->rootJoint()->upperBound (1,  1);
  device->rootJoint()->upperBound (2,  1);

  JointPtr_t ee1 = device->getJointByName ("LAnkleRoll"),
             ee2 = device->getJointByName ("RAnkleRoll"),
             ee3 = device->getJointByName ("RAnklePitch");

  LockedJointPtr_t l1 (new LockedJoint (ee1->rankInConfiguration(), 1, vector_t::Zero(1)));
  LockedJointPtr_t l2 (new LockedJoint (ee2->rankInConfiguration(), 1, vector_t::Zero(1)));
  LockedJointPtr_t l3 (new LockedJoint (ee3->rankInConfiguration(), 1, vector_t::Zero(1)));
  TestFunctionPtr_t t1 (new TestFunction (ee1->rankInConfiguration(), ee2->rankInConfiguration(), 1));
  TestFunctionPtr_t t2 (new TestFunction (ee2->rankInConfiguration(), ee1->rankInConfiguration(), 1));

  RowBlockIndices expectedRow;
  ColBlockIndices expectedCol;

  Configuration_t q = device->currentConfiguration (),
                  qrand = se3::randomConfiguration(device->model());

  {
    ExplicitConstraintSet expression
      (device->configSize(), device->numberDof());
    BOOST_CHECK( expression.add(l1, l1->inArg(), l1->outArg(), l1->inDer(),
                                l1->outDer(), ComparisonTypes_t(1, Equality))
                 >= 0);
    BOOST_CHECK( expression.add(l1, l1->inArg(), l1->outArg(), l1->inDer(),
                                l1->outDer(), ComparisonTypes_t(1, Equality))
                 <  0);
    BOOST_CHECK( expression.add(l2, l2->inArg(), l2->outArg(), l2->inDer(),
                                l2->outDer(), ComparisonTypes_t(1, Equality))
                 >= 0);

    expectedRow = RowBlockIndices();
    expectedRow.addRow (ee1->rankInConfiguration(), 1);
    expectedRow.addRow (ee2->rankInConfiguration(), 1);
    expectedRow.updateRows<true,true,true>();
    BOOST_CHECK_EQUAL (expression.outArgs(), expectedRow);

    expectedRow = RowBlockIndices(BlockIndex::difference (
          BlockIndex::segment_t(0, expression.argSize()),
          expectedRow.rows()));
    BOOST_CHECK_EQUAL (expression.freeArgs(), expectedRow);

    expectedRow = RowBlockIndices();
    expectedRow.addRow (ee1->rankInVelocity(), 1);
    expectedRow.addRow (ee2->rankInVelocity(), 1);
    expectedRow.updateRows<true,true,true>();
    BOOST_CHECK_EQUAL (expression.outDers(), expectedRow);

    expectedCol = ColBlockIndices(BlockIndex::difference (
          BlockIndex::segment_t(0, expression.derSize()),
          expectedRow.rows()));
    BOOST_CHECK_EQUAL (expression.freeDers(), expectedCol);

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
    BOOST_CHECK(expression.viewJacobian (jacobian).eval().isZero());
  }

  {
    ExplicitConstraintSet expression
      (device->configSize(), device->numberDof());
    BOOST_CHECK( expression.add(l1, l1->inArg(), l1->outArg(),
                                l1->inDer(), l1->outDer()) >= 0);
    BOOST_CHECK( expression.add(t1, t1->inArg(), t1->outArg(),
                                t1->inDer(), t1->outDer()) >= 0);

    BOOST_CHECK(expression.solve(qrand));
    vector_t error(expression.outDers().nbIndices());
    BOOST_CHECK(expression.isSatisfied(qrand, error));
    BOOST_CHECK_EQUAL(qrand[ee1->rankInConfiguration()], 0);
    BOOST_CHECK_EQUAL(qrand[ee2->rankInConfiguration()], 0);

    matrix_t jacobian (device->numberDof(), device->numberDof());
    expression.jacobian(jacobian, q);
    BOOST_CHECK(expression.viewJacobian (jacobian).eval().isZero());
  }

  {
    ExplicitConstraintSet expression
      (device->configSize(), device->numberDof());
    BOOST_CHECK( expression.add(t1, t1->inArg(), t1->outArg(), t1->inDer(),
                                t1->outDer()) >= 0);

    matrix_t jacobian (device->numberDof(), device->numberDof());
    expression.jacobian(jacobian, q);
    BOOST_CHECK_EQUAL(jacobian(ee2->rankInVelocity(), ee1->rankInVelocity()), 1);
    BOOST_CHECK_EQUAL(expression.viewJacobian(jacobian).eval().norm(), 1);
  }

  {
    ExplicitConstraintSet expression
      (device->configSize(), device->numberDof());
    BOOST_CHECK( expression.add(t1, t1->inArg(), t1->outArg(), t1->inDer(),
                                t1->outDer()) >= 0);
    BOOST_CHECK( expression.add(t2, t2->inArg(), t2->outArg(), t2->inDer(),
                                t2->outDer()) <  0);
  }

  {
    ExplicitConstraintSet expression
      (device->configSize(), device->numberDof());
    BOOST_CHECK( expression.add(t1, t1->inArg(), t1->outArg(), t1->inDer(),
                                t1->outDer()) >= 0);
    BOOST_CHECK( expression.add(l2, l2->inArg(), l2->outArg(), l2->inDer(),
                                l2->outDer()) <  0);
    BOOST_CHECK( expression.add(l3, l3->inArg(), l3->outArg(), l3->inDer(),
                                l3->outDer()) >= 0);

    matrix_t jacobian (device->numberDof(), device->numberDof());
    expression.jacobian(jacobian, q);
    BOOST_CHECK_EQUAL(jacobian(ee2->rankInVelocity(), ee1->rankInVelocity()), 1);
    BOOST_CHECK_EQUAL(expression.viewJacobian(jacobian).eval().norm(), 1);
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

    ExplicitConstraintSet expression
      (device->configSize(), device->numberDof());
    BOOST_CHECK( expression.add(et, et->inArg(), et->outArg(), et->inDer(),
                                et->outDer()) >= 0);
    BOOST_CHECK( expression.add(l2, l2->inArg(), l2->outArg(), l2->inDer(),
                                l2->outDer()) >= 0);

    matrix_t jacobian (device->numberDof(), device->numberDof());
    expression.jacobian(jacobian, qrand);
  }
}
