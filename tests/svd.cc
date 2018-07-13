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
