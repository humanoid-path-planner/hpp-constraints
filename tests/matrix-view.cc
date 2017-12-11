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

#define EIGEN_RUNTIME_NO_MALLOC

#define BOOST_TEST_MODULE MATRIX_VIEW
#include <boost/test/unit_test.hpp>

#include <iostream>

#include <hpp/constraints/matrix-view.hh>

using namespace Eigen;

// This is an ugly fix to make BOOST_CHECK_EQUAL able to print segments_t
// when they are not equal.
namespace std {
  std::ostream& operator<< (std::ostream& os, BlockIndex::segments_t b)
  {
    internal::print_indices<true>::run (os, b);
    return os;
  }
}

BOOST_AUTO_TEST_CASE(block_index)
{
  typedef BlockIndex::segment_t segment_t;
  typedef BlockIndex::segments_t segments_t;

  segment_t
    a ( 0, 1), // [0]
    b ( 1, 2), // [1,2]
    c ( 0, 0), // []
    d ( 0, 2), // [0,1]
    e ( 4, 3), // [4,6]
    f ( 9, 2), // [9,10]
    g (15, 5); // [15, 19]

  BOOST_CHECK(!BlockIndex::overlap (a, b));
  BOOST_CHECK(!BlockIndex::overlap (a, c));
  BOOST_CHECK(!BlockIndex::overlap (c, b));
  BOOST_CHECK( BlockIndex::overlap (a, a));
  BOOST_CHECK( BlockIndex::overlap (a, d));
  BOOST_CHECK( BlockIndex::overlap (b, d));

  BOOST_CHECK_EQUAL(BlockIndex::difference (a, b), segments_t (1, a));
  BOOST_CHECK_EQUAL(BlockIndex::difference (a, c), segments_t (1, a));
  BOOST_CHECK_EQUAL(BlockIndex::difference (b, d), segments_t (1, segment_t(2,1)));
  BOOST_CHECK_EQUAL(BlockIndex::difference (c, b), segments_t ());
  BOOST_CHECK_EQUAL(BlockIndex::difference (a, a), segments_t ());
  BOOST_CHECK_EQUAL(BlockIndex::difference (a, d), segments_t ());

  segments_t a_plus_f = BlockIndex::sum (a, f);
  BOOST_CHECK_EQUAL(BlockIndex::difference (a_plus_f, b), a_plus_f);

  segments_t v, expected_v, expected_w;
  v.push_back(b);
  v.push_back(a);
  v.push_back(c);
  expected_v.push_back (segment_t(0,3));
  BlockIndex::sort(v);
  BlockIndex::shrink(v);
  BOOST_CHECK_EQUAL(v.size(), 1);
  BOOST_CHECK_EQUAL(BlockIndex::cardinal(v), 3);
  BOOST_CHECK (v == expected_v);

  v.clear();
  BlockIndex::add (v, b);
  BlockIndex::add (v, a);
  BlockIndex::add (v, c);
  BlockIndex::shrink(v);
  BOOST_CHECK_EQUAL(v.size(), 1);
  BOOST_CHECK_EQUAL(BlockIndex::cardinal(v), 3);
  BOOST_CHECK (v == expected_v);

  v.clear ();
  expected_v.clear();
  // v = 0 1 2 3 [4 5 6] 7 8 [9 10] 11 12 13 14 [15 16 17 18 19] 20 ...
  v.push_back (e);
  v.push_back (f);
  v.push_back (g);
  expected_v.push_back (segment_t (5,2));
  expected_v.push_back (f);
  expected_v.push_back (g);
  expected_w.push_back (segment_t (4, 1));
  segments_t w = BlockIndex::split (v, 1);
  BOOST_CHECK (v == expected_v);
  BOOST_CHECK (w == expected_w);

  v.clear ();
  v.push_back (e);
  v.push_back (f);
  v.push_back (g);
  expected_v.clear ();
  expected_v.push_back (segment_t (6,1));
  expected_v.push_back (f);
  expected_v.push_back (g);
  expected_w.clear ();
  expected_w.push_back (segment_t (4, 2));
  w = BlockIndex::split (v, 2);
  BOOST_CHECK (v == expected_v);
  BOOST_CHECK (w == expected_w);

  v.clear ();
  v.push_back (e);
  v.push_back (f);
  v.push_back (g);
  expected_v.clear ();
  expected_v.push_back (f);
  expected_v.push_back (g);
  expected_w.clear ();
  expected_w.push_back (e);
  w = BlockIndex::split (v, 3);
  BOOST_CHECK (v == expected_v);
  BOOST_CHECK (w == expected_w);

  v.clear ();
  v.push_back (e);
  v.push_back (f);
  v.push_back (g);
  expected_v.clear ();
  expected_v.push_back (segment_t (10,1));
  expected_v.push_back (g);
  expected_w.clear ();
  expected_w.push_back (e);
  expected_w.push_back (segment_t (9,1));
  w = BlockIndex::split (v, 4);
  BOOST_CHECK (v == expected_v);
  BOOST_CHECK (w == expected_w);

  v.clear ();
  v.push_back (e);
  v.push_back (f);
  v.push_back (g);
  expected_v.clear ();
  expected_v.push_back (g);
  expected_w.clear ();
  expected_w.push_back (e);
  expected_w.push_back (f);
  w = BlockIndex::split (v, 5);
  BOOST_CHECK (v == expected_v);
  BOOST_CHECK (w == expected_w);

  v.clear ();
  v.push_back (e);
  v.push_back (f);
  v.push_back (g);
  expected_v.clear ();
  expected_w.clear ();
  expected_w.push_back (e);
  expected_w.push_back (f);
  expected_w.push_back (segment_t (15, 1));
  expected_v.push_back (segment_t (16, 4));
  w = BlockIndex::split (v, 6);
  BOOST_CHECK (v == expected_v);
  BOOST_CHECK (w == expected_w);

  v.clear ();
  v.push_back (e);
  v.push_back (f);
  v.push_back (g);
  expected_v.clear ();
  expected_w.clear ();
  expected_w.push_back (e);
  expected_w.push_back (f);
  expected_w.push_back (segment_t (15, 2));
  expected_v.push_back (segment_t (17, 3));
  w = BlockIndex::split (v, 7);
  BOOST_CHECK (v == expected_v);
  BOOST_CHECK (w == expected_w);

  v.clear ();
  v.push_back (e);
  v.push_back (f);
  v.push_back (g);
  expected_v.clear ();
  expected_w.clear ();
  expected_w.push_back (e);
  expected_w.push_back (f);
  expected_w.push_back (segment_t (15, 3));
  expected_v.push_back (segment_t (18, 2));
  w = BlockIndex::split (v, 8);
  BOOST_CHECK (v == expected_v);
  BOOST_CHECK (w == expected_w);

  v.clear ();
  v.push_back (e);
  v.push_back (f);
  v.push_back (g);
  expected_v.clear ();
  expected_w.clear ();
  expected_w.push_back (e);
  expected_w.push_back (f);
  expected_w.push_back (segment_t (15, 4));
  expected_v.push_back (segment_t (19, 1));
  w = BlockIndex::split (v, 9);
  BOOST_CHECK (v == expected_v);
  BOOST_CHECK (w == expected_w);

  v.clear ();
  v.push_back (e);
  v.push_back (f);
  v.push_back (g);
  expected_v.clear ();
  expected_w.clear ();
  expected_w.push_back (e);
  expected_w.push_back (f);
  expected_w.push_back (g);
  w = BlockIndex::split (v, 10);
  BOOST_CHECK (v == expected_v);
  BOOST_CHECK (w == expected_w);

  v.clear ();
  // v = 0 1 2 3 [4 5 6] 7 8 [9 10] 11 12 13 14 [15 16 17 18 19] 20 ...
  v.push_back (e);
  v.push_back (f);
  v.push_back (g);

  expected_w.clear ();
  expected_w.push_back (segment_t (4, 1));
  w = BlockIndex::extract (v, 0, 1);
  BOOST_CHECK (w == expected_w);

  expected_w.clear ();
  expected_w.push_back (segment_t (4, 2));
  w = BlockIndex::extract (v, 0, 2);
  BOOST_CHECK (w == expected_w);

  expected_w.clear ();
  expected_w.push_back (e);
  w = BlockIndex::extract (v, 0, 3);
  BOOST_CHECK (w == expected_w);

  expected_w.clear ();
  expected_w.push_back (e);
  expected_w.push_back (segment_t (9, 1));
  w = BlockIndex::extract (v, 0, 4);
  BOOST_CHECK (w == expected_w);

  expected_w.clear ();
  expected_w.push_back (e);
  expected_w.push_back (f);
  w = BlockIndex::extract (v, 0, 5);
  BOOST_CHECK (w == expected_w);

  expected_w.clear ();
  expected_w.push_back (e);
  expected_w.push_back (f);
  expected_w.push_back (segment_t (15, 1));
  w = BlockIndex::extract (v, 0, 6);
  BOOST_CHECK (w == expected_w);

  expected_w.clear ();
  expected_w.push_back (e);
  expected_w.push_back (f);
  expected_w.push_back (segment_t (15, 2));
  w = BlockIndex::extract (v, 0, 7);
  BOOST_CHECK (w == expected_w);

  expected_w.clear ();
  expected_w.push_back (segment_t (5, 2));
  expected_w.push_back (f);
  expected_w.push_back (segment_t (15, 3));
  w = BlockIndex::extract (v, 1, 7);
  BOOST_CHECK (w == expected_w);

  expected_w.clear ();
  expected_w.push_back (segment_t (6, 1));
  expected_w.push_back (f);
  expected_w.push_back (segment_t (15, 4));
  w = BlockIndex::extract (v, 2, 7);
  BOOST_CHECK (w == expected_w);

  expected_w.clear ();
  expected_w.push_back (f);
  expected_w.push_back (g);
  w = BlockIndex::extract (v, 3, 7);
  BOOST_CHECK (w == expected_w);

  expected_w.clear ();
  expected_w.push_back (f);
  expected_w.push_back (segment_t (15, 4));
  w = BlockIndex::extract (v, 3, 6);
  BOOST_CHECK (w == expected_w);

  expected_w.clear ();
  expected_w.push_back (segment_t (10, 1));
  expected_w.push_back (segment_t (15, 4));
  w = BlockIndex::extract (v, 4, 5);
  BOOST_CHECK (w == expected_w);

  expected_w.clear ();
  expected_w.push_back (segment_t (10, 1));
  expected_w.push_back (segment_t (15, 3));
  w = BlockIndex::extract (v, 4, 4);
  BOOST_CHECK (w == expected_w);
}

template <typename MatrixBlocks_t> void checkMatrixBlocks
(const MatrixBlocks_t& mb, MatrixXd m)
{
  BOOST_CHECK_EQUAL(mb.lview(m).eval(), mb.rview(m).eval());
  BOOST_CHECK_EQUAL(mb.rview(m).eval(),
      mb.rviewTranspose(m.transpose()).eval().transpose());

  MatrixXd res(m);
  res = m;

  mb.lview(res) = mb.rview(m);

  mb.lview(res).setZero();
  BOOST_CHECK (mb.rview(res).isZero());
  BOOST_CHECK (!mb.rview(m).isZero());

  /** CwiseUnaryOp
   *  TODO
   */

  /** CwiseBinaryOp
   *  - check that dynamic allocation is performed.
   */
  MatrixXd res1 (mb.rview(m).rows(), mb.rview(m).cols());
  res = 2 * mb.rview(m).eval();

  Eigen::internal::set_is_malloc_allowed(false);
  // matrix + view
  res1 = mb.rview(m);
  res1 = res1 + mb.rview(m);
  BOOST_CHECK_EQUAL (res1, res);
  // view + matrix
  res1 = mb.rview(m);
  res1 = mb.rview(m) + res1;
  BOOST_CHECK_EQUAL (res1, res);
  // view + view
  MatrixBlocks_t mb2 (mb);
  res1 = mb.rview(m) + mb2.rview(m);
  BOOST_CHECK_EQUAL (res1, res);
  // view of expression
  res1 = mb.rview(m + m);
  BOOST_CHECK_EQUAL (res1, res);
  Eigen::internal::set_is_malloc_allowed(true);

  /** GeneralProduct
   *  TODO
   */
}

BOOST_AUTO_TEST_CASE(matrix_block_view)
{
  typedef MatrixBlocks<false, true> RowsIndices;
  typedef MatrixBlocks<true, false> ColsIndices;
  typedef MatrixBlocks<false, false> MatrixBlocks_t;

  MatrixXd m (10, 11);
  for (MatrixXd::Index i = 0; i < m.rows(); ++i)
    for (MatrixXd::Index j = 0; j < m.cols(); ++j)
      m(i, j) = MatrixXd::Scalar(m.cols() * i + j);

  RowsIndices rows(2,2);
  // rows contains indices 2, 3
  rows.addRow(6, 4);
  // rows contains indices 2, 3, 6, 7, 8, 9

  // Make a ColsIndices from a RowsIndices
  ColsIndices cols (rows);

  MatrixBlocks_t blocks (rows.rows(), cols.cols());

  // Check that operator<< with ostream compiles.
  std::ostringstream oss; oss << rows << '\n' << cols;
  // std::cout << oss.str() << std::endl;

  MatrixXd res, res1;
  rows.lview(m).writeTo(res); // This must resize res.
  BOOST_CHECK_EQUAL(res.rows(), rows.nbRows());
  BOOST_CHECK_EQUAL(res.cols(), m.cols());

  BOOST_CHECK_EQUAL(rows.rview(m).eval().leftCols<8>(), rows.rview(m.leftCols<8>()).eval());

  checkMatrixBlocks (rows, m);
  checkMatrixBlocks (cols, m);
  checkMatrixBlocks (blocks, m);
}
