// Copyright (c) 2018 CNRS
// Authors: Joseph Mirabel
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

#include <hpp/constraints/manipulability.hh>

namespace hpp {
  namespace constraints {
    Manipulability::Manipulability (DifferentiableFunctionPtr_t function,
          DevicePtr_t robot, std::string name) :
      DifferentiableFunction (function->inputSize(),
          function->inputDerivativeSize(), 1, name),
      function_ (function),
      robot_ (robot),
      J_ (function->outputDerivativeSize(), function->inputDerivativeSize())
    {
      activeParameters_           = function->activeParameters();
      activeDerivativeParameters_ = function->activeDerivativeParameters();
      cols_ = Eigen::BlockIndex::fromLogicalExpression (activeDerivativeParameters_);
      J_JT_.resize (J_.rows(),J_.rows());
    }

    void Manipulability::impl_compute (LiegroupElement& res, vectorIn_t arg) const
    {
      assert (cols_.cols().size()>0);

      function_->jacobian (J_, arg);
      value_type logAbsDeterminant;

      // ------------ SVD --------------------------------------------------- //
      J_JT_ = cols_.rview(J_);
      Eigen::JacobiSVD<matrix_t> svd (J_JT_);
      logAbsDeterminant = svd.singularValues().array()
        .cwiseMax(std::numeric_limits<value_type>::min())
        .log10()
        .sum();

      /*
      // ------------ Other decomposition methods --------------------------- //

      // 1. Compute J * J^T
      if (cols_.cols().size() > 1) {
        typedef typename Eigen::ColBlockIndices::View<const matrix_t>::type MatrixView_t;
        MatrixView_t J (cols_.rview(J_));
        //std::cout << J.eval() << std::endl;
        J_JT_.setZero();
        for (MatrixView_t::block_iterator block (J); block.valid(); ++block)
          J_JT_.noalias() += J._block(block) * J._block(block).transpose();
      } else {
        const segment_t& s = cols_.cols()[0];
        //std::cout << J_.middleCols(s.first, s.second) << std::endl;
        J_JT_.noalias() = J_.middleCols(s.first, s.second)
          * J_.middleCols(s.first, s.second).transpose();
      }

      // 2. Compute decomposition

      // 2.1 LDLT
      Eigen::LDLT<matrix_t> ldlt (J_JT_);
      logAbsDeterminant = ldlt.matrixL().nestedExpression().diagonal().array()
        .cwiseMax(std::numeric_limits<value_type>::min())
        .log10()
        .sum();
      logAbsDeterminant += ldlt.vectorD().array()
        .cwiseMax(std::numeric_limits<value_type>::min())
        .log10()
        .sum();
      // 2.2 QRs (FullPiv is more robust that ColPiv)
      //Eigen::ColPivHouseholderQR<matrix_t> qr (J_JT_);
      Eigen::FullPivHouseholderQR<matrix_t> qr (J_JT_);
      logAbsDeterminant = qr.logAbsDeterminant();
      // */

      // This funcion will be used as a cost function whose squared norm is to
      // be minimized.
      res.vector()[0] = std::max(-logAbsDeterminant, 0.);
    }

    void Manipulability::impl_jacobian (matrixOut_t jacobian, vectorIn_t arg) const
    {
      finiteDifferenceCentral (jacobian, arg, robot_, 1e-8);
    }
  } // namespace constraints
} // namespace hpp
