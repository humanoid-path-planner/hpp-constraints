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

#define EIGEN_RUNTIME_NO_MALLOC

#define BOOST_TEST_MODULE MATRIX_VIEW
#include <../tests/util.hh>
#include <boost/test/unit_test.hpp>
#include <hpp/constraints/matrix-view.hh>
#include <iostream>

using namespace Eigen;

BOOST_AUTO_TEST_CASE(block_index) {
  typedef BlockIndex::segment_t segment_t;
  typedef BlockIndex::segments_t segments_t;

  segment_t a(0, 1),  // [0]
      b(1, 2),        // [1,2]
      c(0, 0),        // []
      d(0, 2),        // [0,1]
      e(4, 3),        // [4,6]
      f(9, 2),        // [9,10]
      g(15, 5);       // [15, 19]

  BOOST_CHECK(!BlockIndex::overlap(a, b));
  BOOST_CHECK(!BlockIndex::overlap(a, c));
  BOOST_CHECK(!BlockIndex::overlap(c, b));
  BOOST_CHECK(BlockIndex::overlap(a, a));
  BOOST_CHECK(BlockIndex::overlap(a, d));
  BOOST_CHECK(BlockIndex::overlap(b, d));

  BOOST_CHECK_EQUAL(BlockIndex::difference(a, b), segments_t(1, a));
  BOOST_CHECK_EQUAL(BlockIndex::difference(a, c), segments_t(1, a));
  BOOST_CHECK_EQUAL(BlockIndex::difference(b, d),
                    segments_t(1, segment_t(2, 1)));
  BOOST_CHECK_EQUAL(BlockIndex::difference(c, b), segments_t());
  BOOST_CHECK_EQUAL(BlockIndex::difference(a, a), segments_t());
  BOOST_CHECK_EQUAL(BlockIndex::difference(a, d), segments_t());

  segments_t v, w, expected_v, expected_w;

  v = {a, f};
  BOOST_CHECK_EQUAL(BlockIndex::difference(v, b), v);

  v = {segment_t(0, 5), segment_t(7, 9)};
  v = BlockIndex::difference(v, segment_t(0, 4));
  expected_v = {segment_t(4, 1), segment_t(7, 9)};
  BOOST_CHECK_EQUAL(v, expected_v);

  v = {b, a, c};
  expected_v = {segment_t(0, 3)};
  BlockIndex::sort(v);
  BlockIndex::shrink(v);
  BOOST_CHECK_EQUAL(v.size(), 1);
  BOOST_CHECK_EQUAL(BlockIndex::cardinal(v), 3);
  BOOST_CHECK(v == expected_v);

  v.clear();
  BlockIndex::add(v, b);
  BlockIndex::add(v, a);
  BlockIndex::add(v, c);
  BlockIndex::shrink(v);
  BOOST_CHECK_EQUAL(v.size(), 1);
  BOOST_CHECK_EQUAL(BlockIndex::cardinal(v), 3);
  BOOST_CHECK(v == expected_v);

  w.clear();
  v.clear();
  BlockIndex::add(v, a);
  BlockIndex::add(v, e);
  BlockIndex::add(w, v);

  // v = 0 1 2 3 [4 5 6] 7 8 [9 10] 11 12 13 14 [15 16 17 18 19] 20 ...
  v = {e, f, g};
  expected_v = {segment_t(5, 2), f, g};
  expected_w = {segment_t(4, 1)};
  w = BlockIndex::split(v, 1);
  BOOST_CHECK(v == expected_v);
  BOOST_CHECK(w == expected_w);

  v.clear();
  v = {e, f, g};
  expected_v = {segment_t(6, 1), f, g};
  expected_w = {segment_t(4, 2)};
  w = BlockIndex::split(v, 2);
  BOOST_CHECK(v == expected_v);
  BOOST_CHECK(w == expected_w);

  v = {e, f, g};
  expected_v = {f, g};
  expected_w = {e};
  w = BlockIndex::split(v, 3);
  BOOST_CHECK(v == expected_v);
  BOOST_CHECK(w == expected_w);

  v = {e, f, g};
  expected_v = {segment_t(10, 1), g};
  expected_w = {e, segment_t(9, 1)};
  w = BlockIndex::split(v, 4);
  BOOST_CHECK(v == expected_v);
  BOOST_CHECK(w == expected_w);

  v = {e, f, g};
  expected_v = {g};
  expected_w = {e, f};
  w = BlockIndex::split(v, 5);
  BOOST_CHECK(v == expected_v);
  BOOST_CHECK(w == expected_w);

  v = {e, f, g};
  expected_v = {segment_t(16, 4)};
  expected_w = {e, f, segment_t(15, 1)};
  w = BlockIndex::split(v, 6);
  BOOST_CHECK(v == expected_v);
  BOOST_CHECK(w == expected_w);

  v = {e, f, g};
  expected_v = {segment_t(17, 3)};
  expected_w = {e, f, segment_t(15, 2)};
  w = BlockIndex::split(v, 7);
  BOOST_CHECK(v == expected_v);
  BOOST_CHECK(w == expected_w);

  v = {e, f, g};
  expected_v = {segment_t(18, 2)};
  expected_w = {e, f, segment_t(15, 3)};
  w = BlockIndex::split(v, 8);
  BOOST_CHECK(v == expected_v);
  BOOST_CHECK(w == expected_w);

  v = {e, f, g};
  expected_v = {segment_t(19, 1)};
  expected_w = {e, f, segment_t(15, 4)};
  w = BlockIndex::split(v, 9);
  BOOST_CHECK(v == expected_v);
  BOOST_CHECK(w == expected_w);

  v = {e, f, g};
  expected_v.clear();
  expected_w = {e, f, g};
  w = BlockIndex::split(v, 10);
  BOOST_CHECK(v == expected_v);
  BOOST_CHECK(w == expected_w);

  // v = 0 1 2 3 [4 5 6] 7 8 [9 10] 11 12 13 14 [15 16 17 18 19] 20 ...
  v = {e, f, g};

  expected_w = {segment_t(4, 1)};
  w = BlockIndex::extract(v, 0, 1);
  BOOST_CHECK(w == expected_w);

  expected_w = {segment_t(4, 2)};
  w = BlockIndex::extract(v, 0, 2);
  BOOST_CHECK(w == expected_w);

  expected_w = {e};
  w = BlockIndex::extract(v, 0, 3);
  BOOST_CHECK(w == expected_w);

  expected_w = {e, segment_t(9, 1)};
  w = BlockIndex::extract(v, 0, 4);
  BOOST_CHECK(w == expected_w);

  expected_w = {e, f};
  w = BlockIndex::extract(v, 0, 5);
  BOOST_CHECK(w == expected_w);

  expected_w = {e, f, segment_t(15, 1)};
  w = BlockIndex::extract(v, 0, 6);
  BOOST_CHECK(w == expected_w);

  expected_w = {e, f, segment_t(15, 2)};
  w = BlockIndex::extract(v, 0, 7);
  BOOST_CHECK(w == expected_w);

  expected_w = {segment_t(5, 2), f, segment_t(15, 3)};
  w = BlockIndex::extract(v, 1, 7);
  BOOST_CHECK(w == expected_w);

  expected_w = {segment_t(6, 1), f, segment_t(15, 4)};
  w = BlockIndex::extract(v, 2, 7);
  BOOST_CHECK(w == expected_w);

  expected_w = {f, g};
  w = BlockIndex::extract(v, 3, 7);
  BOOST_CHECK(w == expected_w);

  expected_w = {f, segment_t(15, 4)};
  w = BlockIndex::extract(v, 3, 6);
  BOOST_CHECK(w == expected_w);

  expected_w = {segment_t(10, 1), segment_t(15, 4)};
  w = BlockIndex::extract(v, 4, 5);
  BOOST_CHECK(w == expected_w);

  expected_w = {segment_t(10, 1), segment_t(15, 3)};
  w = BlockIndex::extract(v, 4, 4);
  BOOST_CHECK(w == expected_w);
}

template <typename MatrixBlocks_t>
void checkMatrixBlocks(const MatrixBlocks_t& mb, MatrixXd m) {
  BOOST_TEST_MESSAGE("Current matrix block: " << mb);

  BOOST_CHECK_EQUAL(mb.lview(m).eval(), mb.rview(m).eval());
  BOOST_CHECK_EQUAL(mb.rview(m).eval(),
                    mb.transpose().rview(m.transpose()).eval().transpose());

  MatrixXd res(m);
  res = m;

  mb.lview(res) = mb.rview(m);

  mb.lview(res).setZero();
  BOOST_CHECK(mb.rview(res).isZero());
  BOOST_CHECK(!mb.rview(m).isZero());

  /** Conversion to Ref
   */
#if EIGEN_VERSION_AT_LEAST(3, 2, 92)
  Ref<const MatrixXd> ref(mb.rview(m));
  BOOST_CHECK_EQUAL(ref, mb.rview(m).eval());
  // It is not possible to have a Ref on a MatrixBlockView because
  // Eigen::Ref is based on Eigen::Map.
  // Ref<MatrixXd> ref (mb.lview(m));
#endif  // EIGEN_VERSION_AT_LEAST(3,2,92)

  /** CwiseUnaryOp
   *  TODO
   */

  /** CwiseBinaryOp
   *  - check that dynamic allocation is performed.
   */
  MatrixXd res1(mb.rview(m).rows(), mb.rview(m).cols());
  res = 2 * mb.rview(m).eval();

  Eigen::internal::set_is_malloc_allowed(false);
  // matrix + view
  res1 = mb.rview(m);
  res1 = res1 + mb.rview(m);
  BOOST_CHECK_EQUAL(res1, res);
  // view + matrix
  res1 = mb.rview(m);
  res1 = mb.rview(m) + res1;
  BOOST_CHECK_EQUAL(res1, res);
  // view + view
  MatrixBlocks_t mb2(mb);
  res1 = mb.rview(m) + mb2.rview(m);
  BOOST_CHECK_EQUAL(res1, res);
  // view of expression
  res1 = mb.rview(m + m);
  BOOST_CHECK_EQUAL(res1, res);
  Eigen::internal::set_is_malloc_allowed(true);

  /** GeneralProduct
   *  TODO
   */
}

BOOST_AUTO_TEST_CASE(matrix_block_view) {
  typedef MatrixBlocks<false, true> RowsIndices;
  typedef MatrixBlocks<true, false> ColsIndices;
  typedef MatrixBlocks<false, false> MatrixBlocks_t;

  MatrixXd m(10, 11);
  for (MatrixXd::Index i = 0; i < m.rows(); ++i)
    for (MatrixXd::Index j = 0; j < m.cols(); ++j)
      m(i, j) = MatrixXd::Scalar(m.cols() * i + j);

  RowsIndices rows(2, 2);
  // rows contains indices 2, 3

  // Make a ColsIndices from a RowsIndices
  ColsIndices cols(rows.transpose());

  rows.addRow(6, 4);
  // rows contains indices 2, 3, 6, 7, 8, 9
  cols.addCol(5, 2);
  // cols contains indices 2, 3, 5, 6

  MatrixBlocks_t blocks(rows.rows(), cols.cols());

  // Check that operator<< with ostream compiles.
  std::ostringstream oss;
  oss << rows << '\n' << cols;

  MatrixXd res, res1;
  rows.lview(m).writeTo(res);  // This must resize res.
  BOOST_CHECK_EQUAL(res.rows(), rows.nbRows());
  BOOST_CHECK_EQUAL(res.cols(), m.cols());

  BOOST_CHECK_EQUAL(rows.rview(m).eval().leftCols<8>(),
                    rows.rview(m.leftCols<8>()).eval());

  checkMatrixBlocks(rows, m);
  checkMatrixBlocks(cols, m);
  checkMatrixBlocks(blocks, m);

  checkMatrixBlocks(rows.transpose(), m);
  checkMatrixBlocks(cols.transpose(), m);
  checkMatrixBlocks(blocks.transpose(), m);
  checkMatrixBlocks(blocks.keepRows(), m);
  checkMatrixBlocks(blocks.keepCols(), m);
}

BOOST_AUTO_TEST_CASE(matrix_block_view_iterator) {
  typedef MatrixBlocks<false, true> RowsIndices;
  typedef MatrixBlocks<true, false> ColsIndices;
  typedef MatrixBlocks<false, false> MatrixBlocks_t;

  MatrixXd m(10, 11);
  for (MatrixXd::Index i = 0; i < m.rows(); ++i)
    for (MatrixXd::Index j = 0; j < m.cols(); ++j)
      m(i, j) = MatrixXd::Scalar(m.cols() * i + j);

  RowsIndices rows(2, 2);
  // rows contains indices 2, 3

  // Make a ColsIndices from a RowsIndices
  ColsIndices cols(rows.transpose());

  rows.addRow(6, 4);
  // rows contains indices 2, 3, 6, 7, 8, 9
  cols.addCol(5, 2);
  // cols contains indices 2, 3, 5, 6

  // Extract 2x2 = 4 blocks.
  MatrixBlocks_t blocks(rows.rows(), cols.cols());
  typedef MatrixBlockView<MatrixXd, Dynamic, Dynamic, false, false> MatrixBlockView_t;
  MatrixBlockView_t mbv(blocks.lview(m));
  int i=0;
  typedef MatrixBlockView_t::block_iterator Iterator_t;
  // blocks are traveled by increasing column number. This is quite unexpected.
  for (Iterator_t it(mbv); it.valid(); ++it)
  {
    MatrixXd b;
    switch(i) {
    case 0:
      BOOST_CHECK_EQUAL(it.ri(), 2);
      BOOST_CHECK_EQUAL(it.ci(), 2);
      BOOST_CHECK_EQUAL(it.ro(), 0);
      BOOST_CHECK_EQUAL(it.co(), 0);
      BOOST_CHECK_EQUAL(it.rs(), 2);
      BOOST_CHECK_EQUAL(it.cs(), 2);
      b.resize(2, 2);
      b << 24, 25, 35, 36;
      BOOST_CHECK_EQUAL(mbv._block(it), b);
      break;
    case 1:
      BOOST_CHECK_EQUAL(it.ri(), 6);
      BOOST_CHECK_EQUAL(it.ci(), 2);
      BOOST_CHECK_EQUAL(it.ro(), 2);
      BOOST_CHECK_EQUAL(it.co(), 0);
      BOOST_CHECK_EQUAL(it.rs(), 4);
      BOOST_CHECK_EQUAL(it.cs(), 2);
      b.resize(4, 2);
      b << 68, 69, 79, 80, 90, 91, 101, 102;
      BOOST_CHECK_EQUAL(mbv._block(it), b);
      break;
    case 2:
      BOOST_CHECK_EQUAL(it.ri(), 2);
      BOOST_CHECK_EQUAL(it.ci(), 5);
      BOOST_CHECK_EQUAL(it.ro(), 0);
      BOOST_CHECK_EQUAL(it.co(), 2);
      BOOST_CHECK_EQUAL(it.rs(), 2);
      BOOST_CHECK_EQUAL(it.cs(), 2);
      b.resize(2, 2);
      b << 27, 28, 38, 39;
      BOOST_CHECK_EQUAL(mbv._block(it), b);
      break;
    case 3:
      BOOST_CHECK_EQUAL(it.ri(), 6);
      BOOST_CHECK_EQUAL(it.ci(), 5);
      BOOST_CHECK_EQUAL(it.ro(), 2);
      BOOST_CHECK_EQUAL(it.co(), 2);
      BOOST_CHECK_EQUAL(it.rs(), 4);
      BOOST_CHECK_EQUAL(it.cs(), 2);
      b.resize(4, 2);
      b << 71, 72, 82, 83, 93, 94, 104, 105;
      BOOST_CHECK_EQUAL(mbv._block(it), b);
      break;
    default:
      BOOST_TEST_MESSAGE("This line should not be reached");
      break;
    }
    ++i;
  }

  // Extract 2 rows of full columns
  typedef Eigen::Matrix<double, 10, 11> matrix_10_11_t;
  matrix_10_11_t m1;
  for (auto i = 0; i < m1.rows(); ++i)
    for (auto j = 0; j < m1.cols(); ++j)
      m1(i, j) = MatrixXd::Scalar(m1.cols() * i + j);

  typedef MatrixBlockView<matrix_10_11_t, 10, 11, false, true> MatrixRowView_t;
  MatrixRowView_t mbvCol(rows.lview(m1));
  i=0;
  typedef MatrixRowView_t::block_iterator IteratorRow_t;
  // blocks are traveled by increasing column number. This is quite unexpected.
  for (IteratorRow_t it(mbvCol); it.valid(); ++it)
  {
    MatrixXd b;
    switch(i) {
    case 0:
      BOOST_CHECK_EQUAL(it.ri(), 2);
      BOOST_CHECK_EQUAL(it.ci(), 0);
      BOOST_CHECK_EQUAL(it.ro(), 0);
      BOOST_CHECK_EQUAL(it.co(), 0);
      BOOST_CHECK_EQUAL(it.rs(), 2);
      BOOST_CHECK_EQUAL(it.cs(), 11);
      b.resize(2, 11);
      b << 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
    33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43;
      BOOST_CHECK_EQUAL(mbvCol._block(it), b);
      break;
    case 1:
      BOOST_CHECK_EQUAL(it.ri(), 6);
      BOOST_CHECK_EQUAL(it.ci(), 0);
      BOOST_CHECK_EQUAL(it.ro(), 2);
      BOOST_CHECK_EQUAL(it.co(), 0);
      BOOST_CHECK_EQUAL(it.rs(), 4);
      BOOST_CHECK_EQUAL(it.cs(), 11);
      b.resize(4, 11);
      b << 66,  67,  68,  69,  70,  71,  72,  73,  74, 75,  76,
       77,  78,  79,  80,  81,  82,  83,  84,  85, 86,  87,
       88,  89,  90,  91,  92,  93,  94,  95,  96, 97,  98,
       99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109;
      BOOST_CHECK_EQUAL(mbvCol._block(it), b);
      break;
    default:
      BOOST_TEST_MESSAGE("This line should not be reached");
      break;
    }
    ++i;
  }
  // Extract 2 columns of full rows
  // We use the dynamic size matrix m again in this test.
  typedef MatrixBlockView<MatrixXd, Dynamic, Dynamic, true, false> MatrixColView_t;
  MatrixColView_t mbvRow(cols.lview(m));
  i=0;
  typedef MatrixColView_t::block_iterator IteratorCol_t;
  // blocks are traveled by increasing column number. This is quite unexpected.
  for (IteratorCol_t it(mbvRow); it.valid(); ++it)
  {
    MatrixXd b;
    switch(i) {
    case 0:
      BOOST_CHECK_EQUAL(it.ri(), 0);
      BOOST_CHECK_EQUAL(it.ci(), 2);
      BOOST_CHECK_EQUAL(it.ro(), 0);
      BOOST_CHECK_EQUAL(it.co(), 0);
      BOOST_CHECK_EQUAL(it.rs(), 10);
      BOOST_CHECK_EQUAL(it.cs(), 2);
      b.resize(10, 2);
      b <<  2,   3,
           13,  14,
           24,  25,
           35,  36,
           46,  47,
           57,  58,
           68,  69,
           79,  80,
           90,  91,
          101, 102;
      BOOST_CHECK_EQUAL(mbvRow._block(it), b);
      break;
    case 1:
      BOOST_CHECK_EQUAL(it.ri(), 0);
      BOOST_CHECK_EQUAL(it.ci(), 5);
      BOOST_CHECK_EQUAL(it.ro(), 0);
      BOOST_CHECK_EQUAL(it.co(), 2);
      BOOST_CHECK_EQUAL(it.rs(), 10);
      BOOST_CHECK_EQUAL(it.cs(), 2);
      b.resize(10, 2);
      b <<  5,   6,
           16,  17,
           27,  28,
           38,  39,
           49,  50,
           60,  61,
           71,  72,
           82,  83,
           93,  94,
          104, 105;
      BOOST_CHECK_EQUAL(mbvRow._block(it), b);
      break;
    default:
      BOOST_TEST_MESSAGE("This line should not be reached");
      break;
    }
    ++i;
  }
}
