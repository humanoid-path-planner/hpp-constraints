// Copyright (c) 2018 CNRS
// Authors: Joseph Mirabel
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

    void Manipulability::impl_compute (LiegroupElementRef res, vectorIn_t arg) const
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
