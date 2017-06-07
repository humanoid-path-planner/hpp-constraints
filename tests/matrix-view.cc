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

BOOST_AUTO_TEST_CASE(block_index)
{
  typedef BlockIndex<MatrixXd::Index> BlockIndex_t;
  BlockIndex_t::type
    a ( 0, 1),
    b ( 1, 2),
    c ( 0, 0),
    d ( 0, 2);

  BOOST_CHECK(!BlockIndex_t::overlap (a, b));
  BOOST_CHECK(!BlockIndex_t::overlap (a, c));
  BOOST_CHECK(!BlockIndex_t::overlap (c, b));
  BOOST_CHECK( BlockIndex_t::overlap (a, a));
  BOOST_CHECK( BlockIndex_t::overlap (a, d));
  BOOST_CHECK( BlockIndex_t::overlap (b, d));

  BOOST_CHECK_EQUAL(BlockIndex_t::difference (a, b).size(), 1);
  BOOST_CHECK_EQUAL(BlockIndex_t::difference (a, c).size(), 1);
  BOOST_CHECK_EQUAL(BlockIndex_t::difference (c, b).size(), 0);
  BOOST_CHECK_EQUAL(BlockIndex_t::difference (a, a).size(), 0);
  BOOST_CHECK_EQUAL(BlockIndex_t::difference (a, d).size(), 0);
  BOOST_CHECK_EQUAL(BlockIndex_t::difference (b, d).size(), 1);

  BlockIndex_t::vector_t v;
  v.push_back(b);
  v.push_back(a);
  v.push_back(c);
  BlockIndex_t::sort(v);
  BlockIndex_t::shrink(v);
  BOOST_CHECK_EQUAL(v.size(), 1);
  BOOST_CHECK_EQUAL(BlockIndex_t::cardinal(v), 3);
  BOOST_CHECK(v[0] == BlockIndex_t::type(0, 3));
}

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

BOOST_AUTO_TEST_CASE(matrix_block_view)
{
  typedef MatrixBlockView<const MatrixXd, Dynamic, Dynamic, false, false> MatrixXdConstView;
  typedef MatrixBlockView<MatrixXd, Dynamic, 0, false, true> MatrixRowView;
  typedef MatrixBlockView<VectorXd, Dynamic, 0, false, true> VectorView;

  // EIGEN_STATIC_ASSERT_LVALUE(MatrixRowView)

  typedef MatrixBlockIndexes<false, true> RowsIndexes;
  typedef MatrixBlockIndexes<true, false> ColsIndexes;
  typedef MatrixBlockIndexes<true, false> Indexes;

  MatrixXd m (10, 10);
  for (MatrixXd::Index i = 0; i < m.rows(); ++i)
    for (MatrixXd::Index j = 0; j < m.cols(); ++j)
      m(i, j) = MatrixXd::Scalar(m.cols() * i + j);
  std::cout << m << '\n' << std::endl;

  RowsIndexes rows;
  rows.addRow(2, 2);
  rows.addRow(6, 4);

  ColsIndexes cols;
  cols.addCol(2, 2);
  cols.addCol(6, 4);

  MatrixXd res;
  rows.view(m).writeTo(res);
  std::cout << res << std::endl;

  rows.view(m).setZero();
  std::cout << m << std::endl;

  cols.view(m).writeTo(res);
  std::cout << res << std::endl;

  cols.view(m).setZero();
  std::cout << m << std::endl;

  res = cols.view(m);
  std::cout << res << std::endl;
}
