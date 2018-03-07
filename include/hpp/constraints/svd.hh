// Copyright (c) 2015 CNRS
// Author: Joseph Mirabel
//
//
// This file is part of hpp-constraints
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
// hpp-constraints  If not, see
// <http://www.gnu.org/licenses/>.

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
