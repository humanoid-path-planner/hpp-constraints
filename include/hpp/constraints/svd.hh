// Copyright (c) 2015 CNRS
// Author: Joseph Mirabel
//
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

#ifndef HPP_CONSTRAINTS_SVD_HH
# define HPP_CONSTRAINTS_SVD_HH

# include <hpp/constraints/fwd.hh>
# include <Eigen/SVD>

namespace hpp {
  namespace constraints {
    template <typename SVD>
      static Eigen::Ref<const typename SVD::MatrixUType>
      getU1 (const SVD& svd, const size_type& rank)
    {
      return svd.matrixU().leftCols (rank);
    }

    template <typename SVD>
      static Eigen::Ref<const typename SVD::MatrixUType>
      getU2 (const SVD& svd, const size_type& rank)
    {
      return svd.matrixU().rightCols (svd.matrixU().cols() - rank);
    }

    template <typename SVD>
      static Eigen::Ref<const typename SVD::MatrixUType>
      getV1 (const SVD& svd, const size_type& rank)
    {
      return svd.matrixV().leftCols (rank);
    }

    template <typename SVD>
      static Eigen::Ref<const typename SVD::MatrixUType>
      getV2 (const SVD& svd, const size_type& rank)
    {
      return svd.matrixV().rightCols (svd.matrixV().cols() - rank);
    }

    template < typename SVD>
    static void pseudoInverse(const SVD& svd,
        Eigen::Ref <typename SVD::MatrixType> pinvmat)
    {
      eigen_assert(svd.computeU() && svd.computeV() && "Eigen::JacobiSVD "
          "computation flags must be at least: ComputeThinU | ComputeThinV");

      size_type rank = svd.rank();
      typename SVD::SingularValuesType singularValues_inv =
        svd.singularValues().segment (0,rank).cwiseInverse ();

      pinvmat.noalias() =
        getV1<SVD> (svd, rank) * singularValues_inv.asDiagonal() *
        getU1<SVD> (svd, rank).adjoint();
    }

    template < typename SVD >
    void projectorOnSpan (const SVD& svd,
        Eigen::Ref <typename SVD::MatrixType> projector)
    {
      eigen_assert(svd.computeU() && svd.computeV() && "Eigen::JacobiSVD "
          "computation flags must be at least: ComputeThinU | ComputeThinV");

      size_type rank = svd.rank();
      projector.noalias() = getV1<SVD> (svd, rank) * getV1<SVD>(svd, rank).adjoint();
    }

    template < typename SVD >
    void projectorOnSpanOfInv (const SVD& svd,
        Eigen::Ref <typename SVD::MatrixType> projector)
    {
      eigen_assert(svd.computeU() && svd.computeV() && "Eigen::JacobiSVD "
          "computation flags must be at least: ComputeThinU | ComputeThinV");

      size_type rank = svd.rank();
      projector.noalias() = getU1<SVD>(svd, rank) * getU1<SVD>(svd, rank).adjoint();
    }

    template < typename SVD >
    void projectorOnKernel (const SVD& svd,
        Eigen::Ref <typename SVD::MatrixType> projector,
        const bool& computeFullV = false)
    {
      eigen_assert(svd.computeV() && "Eigen::JacobiSVD "
          "computation flags must be at least: ComputeThinV");

      size_type rank = svd.rank();
      if (computeFullV)
        projector.noalias() = getV2<SVD> (svd, rank) * getV2<SVD>(svd, rank).adjoint();
      else {
        projector.noalias() = - getV1<SVD> (svd, rank) * getV1<SVD>(svd, rank).adjoint();
        projector.diagonal().noalias () += vector_t::Ones(svd.matrixV().rows());
      }
    }

    template < typename SVD >
    void projectorOnKernelOfInv (const SVD& svd,
        Eigen::Ref <typename SVD::MatrixType> projector,
        const bool& computeFullU = false)
    {
      eigen_assert(svd.computeU() && "Eigen::JacobiSVD "
          "computation flags must be at least: ComputeThinU");

      size_type rank = svd.rank();
      if (computeFullU) {
        // U2 * U2*
        projector.noalias() = getU2<SVD>(svd, rank) * getU2<SVD>(svd, rank).adjoint();
      } else {
        // I - U1 * U1*
        projector.noalias() = - getU1<SVD>(svd, rank) * getU1<SVD>(svd, rank).adjoint();
        projector.diagonal().noalias () += vector_t::Ones(svd.matrixU().rows());
      }
    }
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_SVD_HH
