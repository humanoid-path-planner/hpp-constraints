// Copyright (c) 2015, Joseph Mirabel
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
#define BOOST_TEST_MODULE EIGEN_SVD_FEATURE
#include <boost/test/unit_test.hpp>

#include <Eigen/Dense>

#include <hpp/constraints/svd.hh>

using hpp::constraints::pseudoInverse;
using hpp::constraints::projectorOnKernel;
using hpp::constraints::projectorOnKernelOfInv;
using hpp::constraints::projectorOnSpan;
using hpp::constraints::projectorOnSpanOfInv;
using hpp::constraints::matrix_t;
using hpp::constraints::value_type;

template <bool computeFullU, bool computeFullV>
void test ()
{
  const std::size_t rows = 4, cols = 6;
  const value_type tol = 1e-6;
  typedef Eigen::JacobiSVD <matrix_t> SVD;
  int computationFlags =
    ( computeFullU ? Eigen::ComputeFullU : Eigen::ComputeThinU )
    | ( computeFullV ? Eigen::ComputeFullV : Eigen::ComputeThinV );
  SVD svd (rows, cols, computationFlags);
  svd.setThreshold (tol);
  matrix_t Mpinv (cols, rows);
  for (int i = 0; i < 1000; ++i) {
    matrix_t M = matrix_t::Random (rows, cols);
    matrix_t PK (cols, cols);
    matrix_t PS (cols, cols);
    matrix_t PKinv (rows, rows);
    matrix_t PSinv (rows, rows);

    svd.compute (M);
    BOOST_CHECK_MESSAGE ((svd.matrixV().adjoint() * svd.matrixV() - matrix_t::Identity (svd.matrixV().cols(),svd.matrixV().cols())).isZero(),
                         svd.matrixV().adjoint() * svd.matrixV() - matrix_t::Identity (svd.matrixV().cols(),svd.matrixV().cols()));

    pseudoInverse          <SVD> (svd, Mpinv); 

    // TODO There is a multiplication in between two Eigen matrices.
    // I think this is a bug in Eigen.
    // Eigen::internal::set_is_malloc_allowed(false);

    projectorOnKernel      <SVD> (svd, PK, computeFullV); 
    projectorOnSpan        <SVD> (svd, PS);
    projectorOnKernelOfInv <SVD> (svd, PKinv, computeFullU); 
    projectorOnSpanOfInv   <SVD> (svd, PSinv); 

    // Eigen::internal::set_is_malloc_allowed(true);

// This removes a false warning about the conversion sequence used to find the
// proper operator* between matrix_t
#pragma GCC diagnostic ignored "-Wconversion"
    matrix_t Ir = M * Mpinv;
    matrix_t Ic = Mpinv * M;
    matrix_t _M = M * Ic;
    matrix_t _Mpinv = Mpinv * Ir;
#pragma GCC diagnostic pop
    BOOST_CHECK_MESSAGE (_M.isApprox (M), "M = M * M+ * M failed");
    BOOST_CHECK_MESSAGE (_Mpinv.isApprox (Mpinv), "M+ = M+ * M * M+ failed");
    BOOST_CHECK_MESSAGE (Ir.adjoint ().isApprox (Ir), "(M * M+)* = M * M+ failed");
    BOOST_CHECK_MESSAGE (Ic.adjoint ().isApprox (Ic), "(M+ * M)* = M+ * M failed");

    BOOST_CHECK_MESSAGE (PS.isApprox (Ic), "PK = M+ * M failed");
    BOOST_CHECK_MESSAGE (PSinv.isApprox (Ir), "PKinv = M * M+ failed");

    BOOST_CHECK_MESSAGE ((PS + PK)      .isApprox (matrix_t::Identity(cols, cols)), "PS + PK = I failed");
    BOOST_CHECK_MESSAGE ((PSinv + PKinv).isApprox (matrix_t::Identity(rows, rows)), "PSinv + PKinv = I failed");
  }
}

BOOST_AUTO_TEST_CASE(eigen_svd_features)
{
  test<false, false>();
  test<true , false>();
  test<true , true >();
  test<false, true >();
}
