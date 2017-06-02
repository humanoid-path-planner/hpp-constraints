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

#include <hpp/constraints/iterative-solver.hh>

#include <limits>
#include <hpp/util/debug.hh>
#include <hpp/util/timer.hh>

#include <hpp/constraints/svd.hh>
#include <hpp/constraints/macros.hh>

// #define SVD_THRESHOLD Eigen::NumTraits<value_type>::dummy_precision()
#define SVD_THRESHOLD 1e-8

namespace hpp {
  namespace constraints {
    namespace {
      template <bool Superior, bool ComputeJac, typename Derived>
      void compare (value_type& val, const Eigen::MatrixBase<Derived>& jac, const value_type& thr)
      {
        if ((Superior && val < thr) || (!Superior && - thr < val)) {
          if (Superior) val -= thr;
          else          val += thr;
        } else {
          val = 0;
          if (ComputeJac) const_cast<Eigen::MatrixBase<Derived>&> (jac).derived().setZero();
        }
      }

      template <bool ComputeJac>
      void applyComparison (
          const HierarchicalIterativeSolver::ComparisonTypes_t comparison,
          const std::vector<std::size_t>& indexes,
          vector_t& value, matrix_t& jacobian, const value_type& thr)
      {
        for (std::size_t i = 0; i < indexes.size(); ++i) {
          const std::size_t j = indexes[i];
          switch (comparison[j]) {
            case HierarchicalIterativeSolver::Superior: compare<true , ComputeJac> (value[j], jacobian.row(j), thr); break;
            case HierarchicalIterativeSolver::Inferior: compare<false, ComputeJac> (value[j], jacobian.row(j), thr); break;
            default: break;
          }
        }
      }
    }

    HierarchicalIterativeSolver::HierarchicalIterativeSolver()
      : stacks_ (),
      dimension_ (0),
      lastIsOptional_ (false),
      reduction_ (),
      datas_(),
      statistics_ ("HierarchicalIterativeSolver")
    {}

    void HierarchicalIterativeSolver::add (
        const DifferentiableFunctionPtr_t& f,
        const std::size_t& priority,
        const ComparisonTypes_t& comp)
    {
      assert (comp.size() == (std::size_t)f->outputSize());
      const std::size_t minSize = priority + 1;
      if (stacks_.size() < minSize) {
        stacks_.resize (minSize, DifferentiableFunctionStack());
        datas_. resize (minSize, Data());
      }
      stacks_[priority].add(f);
      Data& d = datas_[priority];
      for (std::size_t i = 0; i < comp.size(); ++i) {
        switch (comp[i]) {
          case Superior:
          case Inferior:
            d.inequalityIndexes.push_back (d.comparison.size());
            break;
          case Equality:
            d.equalityIndexes.addRow(d.comparison.size(), 1);
            break;
          default:
            break;
        }
        d.comparison.push_back (comp[i]);
      }
      d.equalityIndexes.updateRows<true, true, true>();
      update();
    }

    void HierarchicalIterativeSolver::update()
    {
      assert(!stacks_.empty());
      // Compute reduced size
      const std::size_t nDofs = stacks_[0].inputDerivativeSize();
      if (reduction_.m_nbCols == 0) reduction_.addCol(0, nDofs);
      std::size_t reducedSize = reduction_.nbIndexes();

      dimension_ = 0;
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        const DifferentiableFunctionStack& f = stacks_[i];
        dimension_ += f.outputSize();
        datas_[i].value        .resize(f.outputSize());
        datas_[i].rightHandSide.resize(f.outputSize());
        datas_[i].rightHandSide.setZero();

        assert((size_type)nDofs == f.inputDerivativeSize());
        datas_[i].jacobian.resize(f.outputDerivativeSize(), f.inputDerivativeSize());
        datas_[i].jacobian.setZero();
        datas_[i].reducedJ.resize(f.outputDerivativeSize(), reducedSize);

        datas_[i].svd = SVD_t (f.outputDerivativeSize(), reducedSize, Eigen::ComputeThinU | Eigen::ComputeThinV);
        datas_[i].svd.setThreshold (SVD_THRESHOLD);
        datas_[i].PK.resize (reducedSize, reducedSize);
      }

      dq_.resize(nDofs);
      dqSmall_.resize(reducedSize);
      projector_.resize(reducedSize, reducedSize);
    }

    vector_t HierarchicalIterativeSolver::rightHandSideFromInput (vectorIn_t arg) const
    {
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        const DifferentiableFunctionStack& f = stacks_[i];
        Data& d = datas_[i];
        f.value (d.value, arg);
        // TODO avoid dynamic allocation
        d.equalityIndexes.view(d.rightHandSide) = d.equalityIndexes.view(d.value).eval();
      }
      return rightHandSide();
    }

    void HierarchicalIterativeSolver::rightHandSide (const vector_t& rhs)
    {
      size_type row = 0;
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        Data& d = datas_[i];
        d.equalityIndexes.view(d.rightHandSide)
          = rhs.segment(row, d.equalityIndexes.m_nbRows);
        row += d.equalityIndexes.m_nbRows;
      }
      assert (row == rhs.size());
    }

    vector_t HierarchicalIterativeSolver::rightHandSide () const
    {
      vector_t rhs;
      size_type row = 0;
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        const Data& d = datas_[i];
        const size_type nRows = d.equalityIndexes.m_nbRows;
        rhs.conservativeResize(rhs.size() + nRows);
        vector_t::SegmentReturnType seg = rhs.segment(row, nRows);
        d.equalityIndexes.view(d.rightHandSide).writeTo(seg);
        row += nRows;
      }
      // assert (row == dimension_);
      return rhs;
    }

    template <bool ComputeJac>
    void HierarchicalIterativeSolver::computeValue (vectorIn_t arg) const
    {
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        const DifferentiableFunctionStack& f = stacks_[i];
        Data& d = datas_[i];

        f.value   (d.value, arg);
        if (ComputeJac) f.jacobian(d.jacobian, arg);
        d.error = d.value - d.rightHandSide;
        applyComparison<ComputeJac>(d.comparison, d.inequalityIndexes, d.error, d.jacobian, inequalityThreshold_);

        // Copy columns that are not reduced
        if (ComputeJac) reduction_.view (d.jacobian).writeTo(d.reducedJ);
      }
    }

    template void HierarchicalIterativeSolver::computeValue<false>(vectorIn_t arg) const;
    template void HierarchicalIterativeSolver::computeValue<true >(vectorIn_t arg) const;

    void HierarchicalIterativeSolver::computeError () const
    {
      const std::size_t end = (lastIsOptional_ ? stacks_.size() - 1 : stacks_.size());
      squaredNorm_ = 0;
      for (std::size_t i = 0; i < end; ++i)
        squaredNorm_ += datas_[i].error.squaredNorm();
    }

    void HierarchicalIterativeSolver::computeDescentDirection () const
    {
      if (stacks_.empty()) {
        dq_.setZero();
        return;
      }
      if (stacks_.size() == 1) { // one level only
        Data& d = datas_[0];
        d.svd.compute (d.reducedJ);
        HPP_DEBUG_SVDCHECK (d.svd);
        // TODO Eigen::JacobiSVD does a dynamic allocation here.
        dqSmall_ = d.svd.solve (- d.error);
      } else {
        projector_.setIdentity();
        vector_t err;
        for (std::size_t i = 0; i < stacks_.size (); ++i) {
          const DifferentiableFunctionStack& f = stacks_[i];
          Data& d = datas_[i];

          // TODO: handle case where this is the first element of the stack and it
          // has no functions
          if (f.outputSize () == 0) continue;
          /// projector is of size numberDof
          bool first = (i == 0);
          bool last = (i == stacks_.size() - 1);
          err = - d.error;
          if (first) {
            // dq should be zero and projector should be identity
            d.svd.compute (d.reducedJ);
            HPP_DEBUG_SVDCHECK (d.svd);
            // TODO Eigen::JacobiSVD does a dynamic allocation here.
            dqSmall_ = d.svd.solve (err);
          } else {
            d.svd.compute (d.reducedJ * projector_);
            HPP_DEBUG_SVDCHECK (d.svd);
            // TODO Eigen::JacobiSVD does a dynamic allocation here.
            dqSmall_ += d.svd.solve (err - d.reducedJ * dqSmall_);
          }
          if (last) break; // No need to compute projector for next step.
          if (!(d.reducedJ * dqSmall_ - err).isZero ()) break;
          /// compute projector for next step.
          projectorOnSpan <SVD_t> (d.svd, d.PK);
          projector_ -= d.PK;
        }
      }
      expandDqSmall();
    }

    void HierarchicalIterativeSolver::expandDqSmall () const
    {
      Eigen::MatrixBlockView<vector_t, Eigen::Dynamic, 1, false, true> (dq_, reduction_.nbIndexes(), reduction_.indexes()) = dqSmall_;
    }
  } // namespace constraints
} // namespace hpp
