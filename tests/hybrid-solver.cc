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
#include <boost/assign/list_of.hpp>

#include <hpp/constraints/hybrid-solver.hh>

#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/simple-device.hh>

#include <hpp/constraints/affine-function.hh>
#include <hpp/constraints/generic-transformation.hh>
#include <hpp/pinocchio/liegroup-element.hh>

#include <../tests/util.hh>

using namespace hpp::constraints;
using boost::assign::list_of;

matrix_t randomPositiveDefiniteMatrix (int N)
{
  matrix_t A (matrix_t::Random(N,N));
  BOOST_REQUIRE ( (A.array() <  1).all() );
  BOOST_REQUIRE ( (A.array() > -1).all() );

  A = (A + A) / 2;
  A += N * matrix_t::Identity (N, N);
  A /= N;
  return A;
}

const value_type test_precision = 1e-5;

//             x       y       z
template <int N1, int N2, int N3>
void test_quadratic ()
{
  const int N = N1 + N2 + N3;

  // (x y z) A (x y z)
  matrix_t A (randomPositiveDefiniteMatrix(N));
  Quadratic::Ptr_t quad (new Quadratic (A));

  // y = B * z
  matrix_t B (matrix_t::Random (N2, N3));
  const int Ninf = std::min(N2,N3);
  B.topLeftCorner (Ninf, Ninf) = randomPositiveDefiniteMatrix(Ninf);
  segment_t in (N1 + N2, N3), out (N1, N2);
  AffineFunctionPtr_t expl (new AffineFunction (B));

  // Make solver
  HybridSolver solver (N, N);
  solver.maxIterations(20);
  solver.errorThreshold(test_precision);
  solver.integration(simpleIntegration<-1,1>);
  solver.saturation(simpleSaturation<-1,1>);

  solver.add (quad, 0);
  solver.explicitConstraintSet().add (expl, in, out, in, out);
  solver.explicitConstraintSetHasChanged();

  matrix_t M (N, N1 + N3);
  M << matrix_t::Identity(N1,N1), matrix_t::Zero(N1,N3),
       matrix_t::Zero    (N2,N1), B,
       matrix_t::Zero    (N3,N1), matrix_t::Identity(N3,N3);
  matrix_t Ar (M.transpose() * A * M);

  BOOST_CHECK_EQUAL (Ar.fullPivLu().rank(), N1 + N3);

  vector_t x (N);

  x.setZero();
  BOOST_CHECK (solver.isSatisfied(x));

  x.setRandom();
  SOLVER_CHECK_SOLVE (solver.solve<solver::lineSearch::Backtracking>(x),
                      SUCCESS);
  // EIGEN_VECTOR_IS_APPROX (x, vector_t::Zero(N));
  EIGEN_VECTOR_IS_APPROX (x.segment<N2>(N1), B * x.tail<N3>());
  BOOST_CHECK_SMALL (value_type(x.transpose() * A * x), test_precision);

  matrix_t expectedJ (1, N1 + N3), J(1, N1 + N3);

  x.setRandom();
  solver.explicitConstraintSet().solve(x);
  expectedJ = 2 * solver.explicitConstraintSet().freeArgs().rview(x).eval().transpose() * Ar;

  solver.computeValue<true> (x);
  solver.updateJacobian(x);
  solver.getReducedJacobian(J);

  EIGEN_IS_APPROX (expectedJ, J);
}

//             w       x       y       z
template <int N1, int N2, int N4, int N3>
void test_quadratic2 ()
{
  const int N = N1 + N2 + N3 + N4;

  // (w x y z) A (w x y z)
  matrix_t A (randomPositiveDefiniteMatrix(N));
  Quadratic::Ptr_t quad (new Quadratic (A));

  // x = B * y
  matrix_t B (matrix_t::Random (N2, N3));
  const int Ninf = std::min(N2,N3);
  B.topLeftCorner (Ninf, Ninf) = randomPositiveDefiniteMatrix(Ninf);
  segment_t in1 (N1 + N2, N3), out1 (N1, N2);
  AffineFunctionPtr_t expl1 (new AffineFunction (B));

  // y = C * z
  matrix_t C (matrix_t::Random (N3, N4));
  const int Ninf2 = std::min(N3,N4);
  C.topLeftCorner (Ninf2, Ninf2) = randomPositiveDefiniteMatrix(Ninf2);
  segment_t in2 (N1 + N2 + N3, N4), out2 (N1 + N2, N3);
  AffineFunctionPtr_t expl2 (new AffineFunction (C));

  // Make solver
  HybridSolver solver (N, N);
  solver.maxIterations(20);
  solver.errorThreshold(test_precision);
  solver.integration(simpleIntegration<-1,1>);
  solver.saturation(simpleSaturation<-1,1>);

  solver.add (quad, 0);
  solver.explicitConstraintSet().add (expl1, in1, out1, in1, out1);
  solver.explicitConstraintSet().add (expl2, in2, out2, in2, out2);
  solver.explicitConstraintSetHasChanged();

  matrix_t M (N, N1 + N4);
  M << matrix_t::Identity(N1,N1), matrix_t::Zero(N1,N4),
       matrix_t::Zero    (N2,N1), B * C,
       matrix_t::Zero    (N3,N1), C,
       matrix_t::Zero    (N4,N1), matrix_t::Identity(N4,N4);
  matrix_t Ar (M.transpose() * A * M);

  BOOST_CHECK_EQUAL (Ar.fullPivLu().rank(), N1 + N4);

  vector_t x (N);

  x.setZero();
  BOOST_CHECK (solver.isSatisfied(x));

  x.setRandom();
  SOLVER_CHECK_SOLVE (solver.solve<solver::lineSearch::Backtracking>(x),
                      SUCCESS);
  // SOLVER_CHECK_SOLVE (solver.solve<lineSearch::Constant>(x), SUCCESS);
  // EIGEN_VECTOR_IS_APPROX (x, vector_t::Zero(N));
  EIGEN_VECTOR_IS_APPROX (x.segment<N2>(N1), B * x.segment<N3>(N1+N2));
  EIGEN_VECTOR_IS_APPROX (x.segment<N3>(N1+N2), C * x.segment<N4>(N1+N2+N3));
  BOOST_CHECK_SMALL (value_type(x.transpose() * A * x), test_precision);

  matrix_t expectedJ (1, N1 + N4), J(1, N1 + N4);

  x.setRandom();
  solver.explicitConstraintSet().solve(x);
  expectedJ = 2 * solver.explicitConstraintSet().freeArgs().rview(x).eval().transpose() * Ar;

  solver.computeValue<true> (x);
  solver.updateJacobian(x);
  solver.getReducedJacobian(J);

  EIGEN_IS_APPROX (expectedJ, J);
}

//             w       x       y       z
template <int N1, int N2, int N4, int N3>
void test_quadratic3 ()
{
  const int N = N1 + N2 + N3 + N4;

  // x = B * (y z)
  matrix_t B (matrix_t::Random (N2, N3 + N4));
  const int Ninf = std::min(N2,N3+N4);
  B.topLeftCorner (Ninf, Ninf) = randomPositiveDefiniteMatrix(Ninf);
  segment_t in1 (N1 + N2, N3 + N4), out1 (N1, N2);
  AffineFunctionPtr_t expl1 (new AffineFunction (B));

  // y = C * z
  matrix_t C (matrix_t::Random (N3, N4));
  const int Ninf2 = std::min(N3,N4);
  C.topLeftCorner (Ninf2, Ninf2) = randomPositiveDefiniteMatrix(Ninf2);
  segment_t in2 (N1 + N2 + N3, N4), out2 (N1 + N2, N3);
  AffineFunctionPtr_t expl2 (new AffineFunction (C));

  // z[0] = d
  vector_t d (vector_t::Random (1));
  segments_t in3; segment_t out3 (N1 + N2 + N3, 1);
  ConstantFunctionPtr_t expl3 (new ConstantFunction (d, 0, 0));

  // (w x y z) A (w x y z)
  matrix_t A (randomPositiveDefiniteMatrix(N));
  Quadratic::Ptr_t quad (new Quadratic (A, -d[0]));

  // Make solver
  HybridSolver solver (N, N);
  solver.maxIterations(20);
  solver.errorThreshold(test_precision);
  solver.integration(simpleIntegration<-1,1>);
  solver.saturation(simpleSaturation<-1,1>);

  solver.add (quad, 0);
  solver.explicitConstraintSet().add (expl1, in1, out1, in1, out1);
  solver.explicitConstraintSet().add (expl2, in2, out2, in2, out2);
  solver.explicitConstraintSet().add (expl3, in3, out3, in3, out3);
  solver.explicitConstraintSetHasChanged();

  matrix_t M (N, N1 + N4);
  M << matrix_t::Identity(N1,N1), matrix_t::Zero(N1,N4),
       matrix_t::Zero    (N2,N1), B.leftCols(N3) * C + B.rightCols(N4),
       matrix_t::Zero    (N3,N1), C,
       matrix_t::Zero    (N4,N1), matrix_t::Identity(N4,N4);
  matrix_t P (N1 + N4, N1 + N4 - 1);
  P << matrix_t::Identity(N1,N1), matrix_t::Zero(N1,N4-1),
       matrix_t::Zero    ( 1,N1+N4-1),
       matrix_t::Zero    (N4-1,N1), matrix_t::Identity(N4-1,N4-1);
  vector_t Xr_0 (vector_t::Zero(N1+N4));
  Xr_0[N1] = d[0];

  matrix_t Ar (M.transpose() * A * M);

  BOOST_CHECK_EQUAL (Ar.fullPivLu().rank(), N1 + N4);

  vector_t x (N);

  x.setRandom();
  SOLVER_CHECK_SOLVE (solver.solve<solver::lineSearch::Backtracking>(x),
                      SUCCESS);
  // SOLVER_CHECK_SOLVE (solver.solve<lineSearch::Constant>(x), SUCCESS);
  // EIGEN_VECTOR_IS_APPROX (x, vector_t::Zero(N));
  EIGEN_VECTOR_IS_APPROX (x.segment<N2>(N1), B * x.segment<N3+N4>(N1+N2));
  EIGEN_VECTOR_IS_APPROX (x.segment<N3>(N1+N2), C * x.segment<N4>(N1+N2+N3));
  BOOST_CHECK_SMALL (value_type(x.transpose() * A * x - d[0]), test_precision);

  matrix_t expectedJ (1, N1 + N4 - 1), J(1, N1 + N4 - 1);

  x.setRandom();
  solver.explicitConstraintSet().solve(x);
  expectedJ = 2 *
    (P * solver.explicitConstraintSet().freeArgs().rview(x).eval() + Xr_0).transpose()
    * Ar * P;

  solver.computeValue<true> (x);
  solver.updateJacobian(x);
  solver.getReducedJacobian(J);

  EIGEN_IS_APPROX (expectedJ, J);
}

BOOST_AUTO_TEST_CASE(quadratic)
{
  test_quadratic<3, 3, 3> ();
  test_quadratic<5, 3, 4> ();

  test_quadratic2<3, 3, 3, 3> ();
  test_quadratic2<3, 4, 2, 6> ();

  test_quadratic3<3, 3, 3, 3> ();
  test_quadratic3<1, 4, 2, 6> ();
}

class LockedJoint : public DifferentiableFunction
{
  public:
    size_type idx_, length_;
    vector_t value_;

    LockedJoint(size_type idx, size_type length, vector_t value)
      : DifferentiableFunction(0, 0, LiegroupSpace::Rn (length), "LockedJoint"),
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

    void impl_compute (LiegroupElement& result, vectorIn_t) const
    {
      result.vector () = value_;
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

    ExplicitTransformation(JointPtr_t joint, size_type in, size_type l,
                           size_type inDer, size_type lDer)
      : DifferentiableFunction(l, lDer,
                               LiegroupSpace::R3xSO3 (),
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
      rt_->value(transform, q);
      result.vector ().head<3>() = transform.vector ().head<3>();
      result.vector ().tail<4>() = Eigen::Quaternion<value_type>
        (exponential(transform.vector ().tail<3>())).coeffs();

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

BOOST_AUTO_TEST_CASE(functions1)
{
  HybridSolver solver(3, 3);

  /// System:
  /// f (q1, q2) = 0
  /// h (    q2) = 0
  ///         q1 = g(q3)
  ///         q2 = C

  // f
  solver.add(AffineFunctionPtr_t(new AffineFunction (matrix_t::Identity(2,3))), 0);
  // q1 = g(q3)
  Eigen::Matrix<value_type,1,1> Jg; Jg (0,0) = 1;
  Eigen::RowBlockIndices inArg; inArg.addRow (2,1);
  Eigen::ColBlockIndices inDer; inDer.addCol (2,1);
  Eigen::RowBlockIndices outArg; outArg.addRow (1,1);
  solver.explicitConstraintSet().add(AffineFunctionPtr_t(new AffineFunction (matrix_t::Ones(1,1))),
      segment_t (2,1), segment_t(0,1),
      segment_t (2,1), segment_t(0,1));
  // q2 = C
  solver.explicitConstraintSet().add(AffineFunctionPtr_t(new AffineFunction (matrix_t(1,0), vector_t::Zero(1))),
      segment_t (), segment_t(1,1),
      segment_t (), segment_t(1,1));

  solver.explicitConstraintSetHasChanged();
  BOOST_CHECK_EQUAL(solver.reducedDimension(), 2);

  // h
  matrix_t h (1,3); h << 0, 1, 0;
  solver.add(AffineFunctionPtr_t(new AffineFunction (h)), 0);
  BOOST_CHECK_EQUAL(solver.       dimension(), 3);
  BOOST_CHECK_EQUAL(solver.reducedDimension(), 2);

  segments_t impDof = list_of(segment_t(2,1));
  BOOST_CHECK_EQUAL(solver.implicitDof(), impDof);
}

BOOST_AUTO_TEST_CASE(functions2)
{
  HybridSolver solver(3, 3);

  /// System:
  /// f (q1, q3) = 0
  /// q2 = g(q3)
  Eigen::Matrix<value_type, 2, 3> Jf;
  Jf << 1, 0, 0,
        0, 0, 1;
  solver.add(AffineFunctionPtr_t(new AffineFunction (Jf)), 0);

  Eigen::Matrix<value_type,1,1> Jg; Jg (0,0) = 1;
  Eigen::RowBlockIndices inArg; inArg.addRow (2,1);
  Eigen::ColBlockIndices inDer; inDer.addCol (2,1);
  Eigen::RowBlockIndices outArg; outArg.addRow (1,1);
  solver.explicitConstraintSet().add(AffineFunctionPtr_t(new AffineFunction (Jg)),
      inArg, outArg, inDer, outArg);

  solver.explicitConstraintSetHasChanged();
  BOOST_CHECK_EQUAL(solver.dimension(), 2);

  // We add to the system h(q3) = 0
  /// f (q1, q3) = 0
  /// h (    q3) = 0
  /// q2 = g(q3)
  // This function should not be removed from the system.
  Eigen::Matrix<value_type, 1, 3> Jh; Jh << 0, 0, 1;
  solver.add(AffineFunctionPtr_t(new AffineFunction (Jh)), 0);
  BOOST_CHECK_EQUAL(solver.dimension(), 3);

  // We add to the system q3 = C
  // Function h should be removed, f should not.
  vector_t C (1); C(0) = 0;
  solver.explicitConstraintSet().add(AffineFunctionPtr_t(new AffineFunction (matrix_t (1, 0), C)),
      segments_t(), segment_t (2, 1),
      segments_t(), segment_t (2, 1));
  solver.explicitConstraintSetHasChanged();

  BOOST_CHECK_EQUAL(solver.       dimension(), 3);
  BOOST_CHECK_EQUAL(solver.reducedDimension(), 2);

  segments_t impDof = list_of(segment_t(0,1));
  BOOST_CHECK_EQUAL(solver.implicitDof(), impDof);
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
  solver.integration(boost::bind(hpp::pinocchio::integrate<true, hpp::pinocchio::DefaultLieGroupMap>, device, _1, _2, _3));
  solver.saturation(boost::bind(saturate, device, _1, _2));

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

  BOOST_CHECK(solver.explicitConstraintSet().add (et, et->inArg(), et->outArg(), et->inDer(), et->outDer()) >= 0);
  solver.explicitConstraintSetHasChanged();
  solver.print(std::cout);

  // BOOST_CHECK_EQUAL(solver.solve<lineSearch::Backtracking  >(q), HybridSolver::SUCCESS);

  Configuration_t tmp = qrand;
  BOOST_CHECK_EQUAL(solver.solve<solver::lineSearch::Backtracking  >(qrand),
                    HybridSolver::SUCCESS);
  qrand = tmp;
  BOOST_CHECK_EQUAL(solver.solve<solver::lineSearch::ErrorNormBased>(qrand),
                    HybridSolver::SUCCESS);
  qrand = tmp;
  BOOST_CHECK_EQUAL(solver.solve<solver::lineSearch::FixedSequence >(qrand),
                    HybridSolver::SUCCESS);

  vector_t dq (device->numberDof());
  dq.setRandom();
  qrand = tmp;
  solver.projectOnKernel (qrand, dq, tmp);
}
