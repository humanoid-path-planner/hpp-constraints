// Copyright (c) 2015, Joseph Mirabel
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

#define BOOST_TEST_MODULE SimplexSolver
#include <boost/test/included/unit_test.hpp>

#include <Eigen/Core>

#include <hpp/constraints/simplex-solver.hh>

using namespace hpp::constraints;
using namespace Eigen;

const Matrix <double, 6, 1> Gravity
= (Matrix <double, 6, 1>() << 0,0,-9.81, 0, 0, 0).finished();

void displayProblem (Ref<MatrixXd> phi, Ref<MatrixXd> A, SimplexSolver& solver) {
  const int n = phi.cols ();
  BOOST_MESSAGE ("\n\n");
  BOOST_MESSAGE ("The problem to solve is phi * F = - m * G, where the unknown "
      "is F >= 0 (a " << n << " dimensional vector).");
  BOOST_MESSAGE ("phi =\n" << phi);

  BOOST_MESSAGE ("Written in LP form: A * (Xp Xm Y)^T = b");
  BOOST_MESSAGE ("A =\n" << A.leftCols (3*n));
  BOOST_MESSAGE ("b^T = " << A.col (3*n).transpose());

  BOOST_MESSAGE ("Optimum: " << solver.getOptimum ());
  BOOST_MESSAGE ("Contact forces : " << (A.col (3*n) - A.leftCols (2*n) * solver.getSolution ().topRows (2*n)).transpose ());
  BOOST_MESSAGE ("             Y = " << solver.getSolution ().segment (2*n, n).transpose());
  BOOST_MESSAGE ("Slack variables: " << solver.getSlack ().transpose());
  BOOST_MESSAGE ("\n\n");
}

BOOST_AUTO_TEST_CASE (buildSignMatrix) {
  Eigen::VectorXd rhs (6);
  Eigen::MatrixXd S (6,6);
  Eigen::MatrixXd res (6,6);
  S.setZero ();
  res.setZero ();
  for (size_t i = 0; i < 100; i++) {
    rhs = Eigen::VectorXd::Random (6);
    SimplexSolver::buildSignMatrix (rhs, S);
    for (int i = 0; i < 6; ++i)
      res(i,i) = (rhs(i) >= 0) ? 1 : -1;
    BOOST_CHECK (S.isApprox (res));
  }
}

BOOST_AUTO_TEST_CASE (solver0) {
  const int n = 2; // Num of contacts
  Eigen::MatrixXd phi    (6,n);
  Eigen::MatrixXd phiinv (n,6);
  phi <<
    0, 0,
    0, 0,
    1, 1,
    0, 0,
    1,-1,
    0, 0;
  phiinv << 0,0, 0.5,0, 0.5,0,
            0,0, 0.5,0,-0.5,0;
  Eigen::MatrixXd I6 (6,6); I6.setIdentity();
  Eigen::MatrixXd In (n,n); In.setIdentity();
  double m = 1;

  Eigen::VectorXd rhs = - m * (phiinv * Gravity);
  Eigen::MatrixXd S (n,n); S.setZero ();
  SimplexSolver::buildSignMatrix (rhs, S);

  Eigen::MatrixXd A (n,n*(2+1)+1);
  A << phiinv * phi - In, In - phiinv * phi, In, rhs;

  SimplexSolver solver (S * A);
  BOOST_CHECK (solver.status () == SimplexSolver::SOLVED);
  displayProblem (phi, A, solver);
}

BOOST_AUTO_TEST_CASE (solver1) {
  const int n = 2; // Num of contacts
  Eigen::MatrixXd phi    (6,n);
  Eigen::MatrixXd phiinv (n,6);
  phi <<
    0,0,
    0,0,
    1,1,
    0,0,
    0.5,2.5,
    0,0;
  phiinv << 0,0, 1.25,0,-0.5,0,
            0,0,-0.25,0, 0.5,0;
  Eigen::MatrixXd I6 (6,6); I6.setIdentity();
  Eigen::MatrixXd In (n,n); In.setIdentity();
  double m = 1;

  Eigen::VectorXd rhs = - m * (phiinv * Gravity);
  Eigen::MatrixXd S (n,n); S.setZero ();
  SimplexSolver::buildSignMatrix (rhs, S);

  Eigen::MatrixXd A (n,n*(2+1)+1);
  A << phiinv * phi - In, In - phiinv * phi, In, rhs;

  SimplexSolver solver (S * A);
  BOOST_CHECK (solver.status () == SimplexSolver::INFEASIBLE);
  displayProblem (phi, A, solver);
}

BOOST_AUTO_TEST_CASE (solver2) {
  const int n = 3; // Num of contacts
  Eigen::MatrixXd phi    (6,n);
  Eigen::MatrixXd phiinv (n,6);
  phi <<
    0, 0, 0,
    0, 0, 0,
    1, 1,-1,
    0, 0, 0,
    1, 3,-2,
    0, 0, 0;
  phiinv << 0,0, 4./3.,0,-0.5,0,
            0,0,-2./3.,0, 0.5,0,
            0,0,-1./3.,0,   0,0;
  Eigen::MatrixXd I6 (6,6); I6.setIdentity();
  Eigen::MatrixXd In (n,n); In.setIdentity();
  double m = 1;

  Eigen::VectorXd rhs = - m * (phiinv * Gravity);
  Eigen::MatrixXd S (n,n); S.setZero ();
  SimplexSolver::buildSignMatrix (rhs, S);

  Eigen::MatrixXd A (n,n*(2+1)+1);
  A << phiinv * phi - In, In - phiinv * phi, In, rhs;

  // The solution are of the form
  // f1 = 2*g +   f2
  // f3 =   g + 2*f2 
  SimplexSolver solver (S * A);
  BOOST_CHECK (solver.status () == SimplexSolver::UNBOUNDED);
  displayProblem (phi, A, solver);
}

BOOST_AUTO_TEST_CASE (solver2bis) {
  const int n = 3; // Num of contacts
  Eigen::MatrixXd phi    (6,n);
  Eigen::MatrixXd phiinv (n,6);
  phi <<
    0, 0, 0,
    0, 0, 0,
    1, 1,-1,
    0, 0, 0,
    1, 3,-2,
    0, 0, 0;
  phiinv << 0,0, 4./3.,0,-0.5,0,
            0,0,-2./3.,0, 0.5,0,
            0,0,-1./3.,0,   0,0;
  Eigen::MatrixXd I6 (6,6); I6.setIdentity();
  Eigen::MatrixXd In (n,n); In.setIdentity();
  double m = 1;

  SimplexSolver solver (3*n, n);
  solver.rhs () = - m * (phiinv * Gravity);
  Eigen::MatrixXd S (n,n); S.setZero ();

  Eigen::MatrixXd A (n,n*(2+1)+1);
  A << phiinv * phi - In, In - phiinv * phi, In, solver.rhs();

  SimplexSolver::buildSignMatrix (solver.rhs(), S);
  solver.rhs().applyOnTheLeft (S);
  solver.constraints () = S * A.leftCols (3*n);

  // The solution are of the form
  // f1 = 2*g +   f2
  // f3 =   g + 2*f2 
  solver.checkFeasability ();
  solver.buildSolution ();

  BOOST_CHECK (solver.status () == SimplexSolver::UNBOUNDED);
  displayProblem (phi, A, solver);
}

BOOST_AUTO_TEST_CASE (solver3) {
  const int n = 4; // Num of contacts
  Eigen::MatrixXd phi    (6,n);
  Eigen::MatrixXd phiinv (n,6);
  phi <<
    0, 0, 0, 0,
    0, 0, 0, 0,
    1, 1,-1,-1,
   -1,-2, 1, 2,
    1, 1,-1,-1,
    0, 0, 0, 0;
  //phiinv << 0,0, 0.75,0,-0.25,0,
            //0,0,-0.25,0, 0.25,0,
            //0,0,-0.75,0, 0.25,0,
            //0,0, 0.25,0,-0.25,0;
  phiinv << 0,0, 0.5 , 0.5, 0.5 ,0,
            0,0,-0.25,-0.5,-0.25,0,
            0,0,-0.5 ,-0.5,-0.5 ,0,
            0,0, 0.25, 0.5, 0.25,0;
  Eigen::MatrixXd I6 (6,6); I6.setIdentity();
  Eigen::MatrixXd In (n,n); In.setIdentity();
  double m = 1;

  Eigen::VectorXd rhs = - m * (phiinv * Gravity);
  Eigen::MatrixXd S (n,n); S.setZero ();
  SimplexSolver::buildSignMatrix (rhs, S);

  Eigen::MatrixXd A (n,n*(2+1)+1);
  A << phiinv * phi - In, In - phiinv * phi, In, rhs;

  // The solution are of the form
  SimplexSolver solver (S * A);
  BOOST_CHECK (solver.status () == SimplexSolver::UNBOUNDED);
  displayProblem (phi, A, solver);
}
