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

#define BOOST_TEST_MODULE MATRIX_VIEW
#include <boost/test/unit_test.hpp>

#include <iostream>

#include <hpp/constraints/matrix-view.hh>

using namespace Eigen;

BOOST_AUTO_TEST_CASE(matrix_view)
{
  typedef MatrixView<const MatrixXd, Dynamic, Dynamic, false, false> MatrixXdConstView;
  typedef MatrixView<MatrixXd, Dynamic, 0, false, true> MatrixRowView;
  typedef MatrixView<VectorXd, Dynamic, 0, false, true> VectorView;

  EIGEN_STATIC_ASSERT_LVALUE(MatrixRowView)

  MatrixXd m (10, 10);
  for (MatrixXd::Index i = 0; i < m.rows(); ++i)
    for (MatrixXd::Index j = 0; j < m.cols(); ++j)
      m(i, j) = MatrixXd::Scalar(m.cols() * i + j);
  std::cout << m << '\n' << std::endl;

  MatrixRowView::Indexes_t rows(5);
  rows[0] = 0;
  rows[1] = 2;
  rows[2] = 3;
  rows[3] = 5;
  rows[4] = 9;

  MatrixXd ret = MatrixRowView(m, rows);
  std::cout << ret << '\n' << std::endl;

  VectorXd x (VectorXd::Random(m.cols()));
  VectorXd y1 = m   * x;
  VectorXd y2 = ret * x;
  VectorXd y3 = MatrixRowView(m, rows) * x;

  BOOST_CHECK(y3.isApprox(y2));
  for (std::size_t i = 0; i < rows.size(); ++i) {
    BOOST_CHECK_CLOSE(y1(rows[i]), y2(i), NumTraits<double>::dummy_precision());
  }
  BOOST_CHECK (VectorView (y1, rows).isApprox(y2));

  // Read-only access
  MatrixXdConstView::Indexes_t cols(10);
  for (std::size_t i = 0; i < cols.size(); ++i) cols[i] = 9 - i;
  std::cout << MatrixXdConstView(m, rows, cols) << '\n' << std::endl;

  // Write access
  MatrixRowView(m, rows).setZero();
  std::cout << m << '\n' << std::endl;
}

