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
#include <hpp/constraints/impl/iterative-solver.hh>

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

    namespace lineSearch {
      template bool Constant::operator() (const HierarchicalIterativeSolver& solver, vectorOut_t arg, vectorOut_t darg);

      Backtracking::Backtracking () : c (0.001), tau (0.7), smallAlpha (0.2) {}
      template bool Backtracking::operator() (const HierarchicalIterativeSolver& solver, vectorOut_t arg, vectorOut_t darg);

      FixedSequence::FixedSequence() : alpha (.2), alphaMax (.95), K (.8) {}
      template bool FixedSequence::operator() (const HierarchicalIterativeSolver& solver, vectorOut_t arg, vectorOut_t darg);

      ErrorNormBased::ErrorNormBased(value_type alphaMin, value_type _a, value_type _b)
          : C (0.5 + alphaMin / 2), K ((1 - alphaMin) / 2), a (_a), b (_b)
      {}
      template bool ErrorNormBased::operator() (const HierarchicalIterativeSolver& solver, vectorOut_t arg, vectorOut_t darg);
    }

    HierarchicalIterativeSolver::HierarchicalIterativeSolver (const std::size_t& argSize, const std::size_t derSize)
      : stacks_ (),
      argSize_ (argSize),
      derSize_ (derSize),
      dimension_ (0),
      lastIsOptional_ (false),
      reduction_ (),
      datas_(),
      statistics_ ("HierarchicalIterativeSolver")
    {
      reduction_.addCol (0, derSize_);
    }

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

    bool_array_t HierarchicalIterativeSolver::activeParameters () const
    {
      bool_array_t ap (bool_array_t::Constant(argSize_, false));
      for (std::size_t i = 0; i < stacks_.size (); ++i)
        ap = ap || stacks_[i].activeParameters();
      return ap;
    }

    bool_array_t HierarchicalIterativeSolver::activeDerivativeParameters () const
    {
      bool_array_t ap (bool_array_t::Constant(derSize_, false));
      for (std::size_t i = 0; i < stacks_.size (); ++i)
        ap = ap || stacks_[i].activeDerivativeParameters();
      return ap;
    }

    void HierarchicalIterativeSolver::update()
    {
      // Compute reduced size
      std::size_t reducedSize = reduction_.nbIndexes();

      dimension_ = 0;
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        computeActiveRowsOfJ (i);

        const DifferentiableFunctionStack& f = stacks_[i];
        dimension_ += datas_[i].activeRowsOfJ.m_nbRows;
        datas_[i].output = LiegroupElement (f.outputSpace ());
        datas_[i].rightHandSide = LiegroupElement (f.outputSpace ());
        datas_[i].rightHandSide.setNeutral ();

        assert(derSize_ == f.inputDerivativeSize());
        datas_[i].jacobian.resize(f.outputDerivativeSize(), f.inputDerivativeSize());
        datas_[i].jacobian.setZero();
        datas_[i].reducedJ.resize(datas_[i].activeRowsOfJ.m_nbRows, reducedSize);

        datas_[i].svd = SVD_t (f.outputDerivativeSize(), reducedSize, Eigen::ComputeThinU | Eigen::ComputeThinV);
        datas_[i].svd.setThreshold (SVD_THRESHOLD);
        datas_[i].PK.resize (reducedSize, reducedSize);

        datas_[i].maxRank = 0;
      }

      dq_ = vector_t::Zero(derSize_);
      dqSmall_.resize(reducedSize);
      projector_.resize(reducedSize, reducedSize);
      reducedJ_.resize(dimension_, reducedSize);
      svd_ = SVD_t (dimension_, reducedSize, Eigen::ComputeThinU | Eigen::ComputeThinV);
    }

    void HierarchicalIterativeSolver::computeActiveRowsOfJ (std::size_t iStack)
    {
      Data& d = datas_[iStack];
      const DifferentiableFunctionStack& f = stacks_[iStack];
      const DifferentiableFunctionStack::Functions_t& fs = f.functions();
      std::size_t row = 0;

      typedef Eigen::MatrixBlocks<false, false> BlockIndexes;
      BlockIndexes::segments_t rows;
      for (std::size_t i = 0; i < fs.size (); ++i) {
        bool_array_t adp = reduction_.rviewTranspose(fs[i]->activeDerivativeParameters().matrix()).eval();
        if (adp.any()) // If at least one element of adp is true
          rows.push_back (BlockIndexes::segment_t
                          (row, fs[i]->outputDerivativeSize()));
        row += fs[i]->outputDerivativeSize();
      }
      d.activeRowsOfJ = Eigen::MatrixBlocks<false,false> (rows, reduction_.m_cols);
      d.activeRowsOfJ.updateRows<true, true, true>();
    }

    vector_t HierarchicalIterativeSolver::rightHandSideFromInput (vectorIn_t arg) const
    {
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        const DifferentiableFunctionStack& f = stacks_[i];
        Data& d = datas_[i];
        f.value (d.output, arg);
        // TODO avoid dynamic allocation
        d.equalityIndexes.lview(d.rightHandSide.vector ()) =
          d.equalityIndexes.rview(d.output.vector ()).eval();
      }
      return rightHandSide();
    }

    bool HierarchicalIterativeSolver::rightHandSideFromInput (const DifferentiableFunctionPtr_t& f, vectorIn_t arg) const
    {
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        Data& d = datas_[i];
        const DifferentiableFunctionStack::Functions_t& fs = stacks_[i].functions();
        size_type row = 0;
        for (std::size_t j = 0; j < fs.size(); ++j) {
          if (f == fs[j]) {
            LiegroupElement tmp (f->outputSpace ());
            f->value (tmp, arg);
            d.output.vector ().segment(row, f->outputSize()) = tmp.vector ();
            for (size_type k = 0; k < f->outputSize(); ++k) {
              if (d.comparison[row + k] == Equality) {
                d.rightHandSide.vector () [row + k] = d.output.vector ()[row + k];
              }
            }
            return true;
          }
          row += fs[i]->outputSize();
        }
      }
      return false;
    }

    bool HierarchicalIterativeSolver::rightHandSide (const DifferentiableFunctionPtr_t& f, vectorIn_t rhs) const
    {
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        Data& d = datas_[i];
        const DifferentiableFunctionStack::Functions_t& fs = stacks_[i].functions();
        size_type row = 0;
        for (std::size_t j = 0; j < fs.size(); ++j) {
          if (f == fs[j]) {
            for (size_type k = 0; k < f->outputSize(); ++k) {
              if (d.comparison[row + k] == Equality) {
                d.rightHandSide.vector () [row + k] = rhs [row + k];
              }
            }
            return true;
          }
          row += fs[i]->outputSize();
        }
      }
      return false;
    }

    void HierarchicalIterativeSolver::rightHandSide (vectorIn_t rhs)
    {
      size_type row = 0;
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        Data& d = datas_[i];
        d.equalityIndexes.lview(d.rightHandSide.vector ())
          = rhs.segment(row, d.equalityIndexes.m_nbRows);
        row += d.equalityIndexes.m_nbRows;
      }
      assert (row == rhs.size());
    }

    vector_t HierarchicalIterativeSolver::rightHandSide () const
    {
      vector_t rhs(rightHandSideSize());
      size_type row = 0;
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        const Data& d = datas_[i];
        const size_type nRows = d.equalityIndexes.m_nbRows;
        vector_t::SegmentReturnType seg = rhs.segment(row, nRows);
        seg = d.equalityIndexes.rview(d.rightHandSide.vector ());
        row += nRows;
      }
      assert (row == rhs.size());
      return rhs;
    }

    size_type HierarchicalIterativeSolver::rightHandSideSize () const
    {
      size_type rhsSize = 0;
      for (std::size_t i = 0; i < stacks_.size (); ++i)
        rhsSize += datas_[i].equalityIndexes.m_nbRows;
      return rhsSize;
    }

    template <bool ComputeJac>
    void HierarchicalIterativeSolver::computeValue (vectorIn_t arg) const
    {
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        const DifferentiableFunctionStack& f = stacks_[i];
        Data& d = datas_[i];

        f.value   (d.output, arg);
        if (ComputeJac) f.jacobian(d.jacobian, arg);
        d.error = d.output - d.rightHandSide;
        applyComparison<ComputeJac>(d.comparison, d.inequalityIndexes, d.error, d.jacobian, inequalityThreshold_);

        // Copy columns that are not reduced
        if (ComputeJac) d.reducedJ = d.activeRowsOfJ.rview (d.jacobian);
      }
    }

    template void HierarchicalIterativeSolver::computeValue<false>(vectorIn_t arg) const;
    template void HierarchicalIterativeSolver::computeValue<true >(vectorIn_t arg) const;

    void HierarchicalIterativeSolver::getValue (vectorOut_t v) const
    {
      size_type row = 0;
      for (std::size_t i = 0; i < datas_.size(); ++i) {
        const Data& d = datas_[i];
        v.segment(row, d.output.vector ().rows()) = d.output.vector ();
        row += d.output.vector ().rows();
      }
    }

    void HierarchicalIterativeSolver::getReducedJacobian (matrixOut_t J) const
    {
      size_type row = 0;
      for (std::size_t i = 0; i < datas_.size(); ++i) {
        const Data& d = datas_[i];
        J.middleRows(row, d.reducedJ.rows()) = d.reducedJ;
        row += d.reducedJ.rows();
      }
    }

    void HierarchicalIterativeSolver::computeError () const
    {
      const std::size_t end = (lastIsOptional_ ? stacks_.size() - 1 : stacks_.size());
      squaredNorm_ = 0;
      for (std::size_t i = 0; i < end; ++i)
        squaredNorm_ += datas_[i].error.squaredNorm();
    }

    void HierarchicalIterativeSolver::computeDescentDirection () const
    {
      sigma_ = std::numeric_limits<value_type>::max();

      if (stacks_.empty()) {
        dq_.setZero();
        return;
      }
      if (stacks_.size() == 1) { // one level only
        Data& d = datas_[0];
        d.svd.compute (d.reducedJ);
        HPP_DEBUG_SVDCHECK (d.svd);
        // TODO Eigen::JacobiSVD does a dynamic allocation here.
        dqSmall_ = d.svd.solve (- Eigen::RowBlockIndexes(d.activeRowsOfJ.m_rows).rview(d.error).eval());
        d.maxRank = std::max(d.maxRank, d.svd.rank());
        if (d.maxRank > 0)
          sigma_ = std::min(sigma_, d.svd.singularValues()[d.maxRank - 1]);
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
          err = - Eigen::RowBlockIndexes(d.activeRowsOfJ.m_rows).rview(d.error).eval();
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
          // Update sigma
          d.maxRank = std::max(d.maxRank, d.svd.rank());
          if (d.maxRank > 0)
            sigma_ = std::min(sigma_, d.svd.singularValues()[d.maxRank - 1]);

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

    template HierarchicalIterativeSolver::Status HierarchicalIterativeSolver::solve (vectorOut_t arg, lineSearch::Backtracking   lineSearch) const;
    template HierarchicalIterativeSolver::Status HierarchicalIterativeSolver::solve (vectorOut_t arg, lineSearch::FixedSequence  lineSearch) const;
    template HierarchicalIterativeSolver::Status HierarchicalIterativeSolver::solve (vectorOut_t arg, lineSearch::ErrorNormBased lineSearch) const;
  } // namespace constraints
} // namespace hpp
