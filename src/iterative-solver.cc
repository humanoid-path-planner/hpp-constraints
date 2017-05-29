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
      HPP_DEFINE_TIMECOUNTER (iterative_solver);
    }

    HPP_DEFINE_REASON_FAILURE (REASON_MAX_ITER, "Max Iterations reached");
    HPP_DEFINE_REASON_FAILURE (REASON_ERROR_INCREASED, "Error increased");

    HierarchicalIterativeSolver::HierarchicalIterativeSolver()
      : Solver (Solver::Iterative),
      stacks_ (),
      dimension_ (0),
      lastIsOptional_ (false),
      datas_(),
      statistics_ ("HierarchicalIterativeSolver")
    {}

    void HierarchicalIterativeSolver::update()
    {
      assert(!stacks_.empty());
      // Compute reduced size
      const std::size_t nDofs = stacks_[0].inputDerivativeSize();
      std::size_t reducedSize = 0;
      if (reduction_.empty()) reducedSize = nDofs;
      for (intervals_t::const_iterator _int = reduction_.begin ();
          _int != reduction_.end (); ++_int)
        reducedSize += _int->second;

      dimension_ = 0;
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        const DifferentiableFunctionStack& f = stacks_[i];
        dimension_ += f.outputSize();
        datas_[i].value        .resize(f.outputSize());
        datas_[i].rightHandSide.resize(f.outputSize());
        datas_[i].rightHandSide.setZero();

        assert(nDofs == f.inputDerivativeSize());
        datas_[i].jacobian.resize(f.outputDerivativeSize(), f.inputDerivativeSize());
        datas_[i].jacobian.setZero();
        datas_[i].reducedJ.resize(f.outputDerivativeSize(), reducedSize);

        datas_[i].svd = SVD_t (f.outputSize(), reducedSize, Eigen::ComputeThinU | Eigen::ComputeThinV);
        datas_[i].svd.setThreshold (SVD_THRESHOLD);
        datas_[i].PK.resize (reducedSize, reducedSize);
      }

      dq_.resize(nDofs);
      dqSmall_.resize(reducedSize);
      projector_.resize(reducedSize, reducedSize);
    }

    bool HierarchicalIterativeSolver::solve (vectorOut_t arg) const
    {
      // hppDout (info, "before projection: " << arg.transpose ());
      assert (!arg.hasNaN());
      HPP_START_TIMECOUNTER (iterative_solver);

      value_type alpha = .2;
      value_type alphaMax = .95;
      size_type errorDecreased = 3, iter = 0;
      value_type previousSquaredNorm =
	std::numeric_limits<value_type>::infinity();

      // Fill value and Jacobian
      computeValueAndJacobian (arg);
      computeError();

      while (squaredNorm_ > squaredErrorThreshold_ && errorDecreased &&
	     iter < maxIterations_) {

        computeIncrement (alpha);
        integrate_(arg, dq_, arg);

	// Increase alpha towards alphaMax
	computeValueAndJacobian (arg);
	alpha = alphaMax - .8*(alphaMax - alpha);
        computeError ();
	hppDout (info, "squareNorm = " << squaredNorm_);
	--errorDecreased;
	if (squaredNorm_ < previousSquaredNorm) errorDecreased = 3;
	previousSquaredNorm = squaredNorm_;
	++iter;

      }

      if (squaredNorm_ > squaredErrorThreshold_) {
        statistics_.addFailure ((!errorDecreased)?REASON_ERROR_INCREASED:REASON_MAX_ITER);
        statistics_.isLowRatio (true);
      } else {
        statistics_.addSuccess();
      }
      HPP_STOP_TIMECOUNTER (iterative_solver);
      HPP_DISPLAY_TIMECOUNTER (iterative_solver);
      hppDout (info, "number of iterations: " << iter);
      if (squaredNorm_ > squaredErrorThreshold_) {
	hppDout (info, "Projection failed.");
	return false;
      }
      hppDout (info, "After projection: " << arg.transpose ());
      assert (!arg.hasNaN());
      return true;
    }

    vector_t HierarchicalIterativeSolver::rightHandSideFromInput (vectorIn_t arg) const
    {
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        const DifferentiableFunctionStack& f = stacks_[i];
        Data& d = datas_[i];
        f.value (d.rightHandSide, arg);
      }
      return rightHandSide();
    }

    void HierarchicalIterativeSolver::rightHandSide (const vector_t& rhs)
    {
      assert (rhs.size() == dimension_);
      size_type row = 0;
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        const DifferentiableFunctionStack& f = stacks_[i];
        Data& d = datas_[i];
        d.rightHandSide = rhs.segment (row, d.rightHandSide.size());
        row += d.rightHandSide.size();
      }
      assert (row == dimension_);
    }

    vector_t HierarchicalIterativeSolver::rightHandSide () const
    {
      vector_t rhs(dimension_);
      size_type row = 0;
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        const DifferentiableFunctionStack& f = stacks_[i];
        const Data& d = datas_[i];
        rhs.segment(row, d.rightHandSide.size()) = d.rightHandSide;
        row += d.rightHandSide.size();
      }
      assert (row == dimension_);
      return rhs;
    }

    void HierarchicalIterativeSolver::computeValueAndJacobian (vectorIn_t arg, const std::size_t priority) const
    {
      assert(priority < stacks_.size());
      const DifferentiableFunctionStack& f = stacks_[priority];
      Data& d = datas_[priority];

      f.value   (d.value, arg);
      f.jacobian(d.jacobian, arg);
      // TODO (*(*it)->comparisonType ()) (v, jacobian);
      d.error = d.value - d.rightHandSide;

      size_type col = 0;
      // Copy columns that are not reduced
      for (intervals_t::const_iterator _int = reduction_.begin ();
          _int != reduction_.end (); ++_int) {
        size_type col0 = _int->first;
        size_type nbCols = _int->second;
        d.reducedJ.middleCols (col, nbCols) =
          d.jacobian.middleCols (col0, nbCols);
        col += nbCols;
      }
    }

    inline void HierarchicalIterativeSolver::computeValueAndJacobian (vectorIn_t arg) const
    {
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        computeValueAndJacobian(arg, i);
      }
    }

    inline void HierarchicalIterativeSolver::computeError () const
    {
      const std::size_t end = (lastIsOptional_ ? stacks_.size() - 1 : stacks_.size());
      squaredNorm_ = 0;
      for (std::size_t i = 0; i < end; ++i) {
        datas_[i].error = datas_[i].value - datas_[i].rightHandSide;
        squaredNorm_ += datas_[i].value.squaredNorm();
      }
    }

    void HierarchicalIterativeSolver::computeIncrement (const value_type& alpha) const
    {
      if (stacks_.empty()) {
        dq_.setZero();
        return;
      }
      if (stacks_.size() == 1) { // one level only
        Data& d = datas_[0];
        d.svd.compute (d.reducedJ);
        HPP_DEBUG_SVDCHECK (d.svd);
        dqSmall_ = d.svd.solve (- alpha * d.error);
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
          err = - alpha * d.error;
          if (first) {
            // dq should be zero and projector should be identity
            d.svd.compute (d.reducedJ);
            HPP_DEBUG_SVDCHECK (d.svd);
            dqSmall_ = d.svd.solve (err);
          } else {
            d.svd.compute (d.reducedJ * projector_);
            HPP_DEBUG_SVDCHECK (d.svd);
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
      if (reduction_.empty ()) {
	dq_ = dqSmall_;
	return;
      }
      size_type col = 0;
      for (intervals_t::const_iterator _int = reduction_.begin ();
          _int != reduction_.end (); ++_int) {
	size_type col0 = _int->first;
	size_type nbCols = _int->second;

	dq_.segment (col0, nbCols) = dqSmall_.segment (col, nbCols);
	col += nbCols;
      }
    }
  } // namespace constraints
} // namespace hpp
