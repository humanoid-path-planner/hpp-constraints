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

#include <hpp/constraints/solver/hierarchical-iterative.hh>
#include <hpp/constraints/solver/impl/hierarchical-iterative.hh>

#include <limits>
#include <hpp/util/debug.hh>
#include <hpp/util/timer.hh>

#include <hpp/pinocchio/util.hh>
#include <hpp/pinocchio/liegroup-element.hh>

#include <hpp/constraints/svd.hh>
#include <hpp/constraints/macros.hh>
#include <hpp/constraints/implicit.hh>

// #define SVD_THRESHOLD Eigen::NumTraits<value_type>::dummy_precision()
#define SVD_THRESHOLD 1e-8

namespace hpp {
  namespace constraints {
    namespace solver {
      namespace {
        template <bool Superior, bool ComputeJac, typename Derived>
        void compare (value_type& val, const Eigen::MatrixBase<Derived>& jac,
                      const value_type& thr)
        {
          if ((Superior && val < thr) || (!Superior && - thr < val)) {
            if (Superior) val -= thr;
            else          val += thr;
          } else {
            val = 0;
            if (ComputeJac) const_cast<Eigen::MatrixBase<Derived>&>
                              (jac).derived().setZero();
          }
        }

        template <bool ComputeJac>
        void applyComparison (const ComparisonTypes_t comparison,
                              const std::vector<std::size_t>& indices,
                              vector_t& value, matrix_t& jacobian,
                              const value_type& thr)
        {
          for (std::size_t i = 0; i < indices.size(); ++i) {
            const std::size_t j = indices[i];
            switch (comparison[j]) {
            case Superior: compare<true , ComputeJac>
                (value[j], jacobian.row(j), thr); break;
            case Inferior: compare<false, ComputeJac>
                (value[j], jacobian.row(j), thr); break;
            default: break;
            }
          }
        }
      }

      namespace lineSearch {
        template bool Constant::operator()
          (const HierarchicalIterative& solver, vectorOut_t arg, vectorOut_t darg);

        Backtracking::Backtracking () : c (0.001), tau (0.7), smallAlpha (0.2) {}
        template bool Backtracking::operator()
          (const HierarchicalIterative& solver, vectorOut_t arg, vectorOut_t darg);

        FixedSequence::FixedSequence() : alpha (.2), alphaMax (.95), K (.8) {}
        template bool FixedSequence::operator()
          (const HierarchicalIterative& solver, vectorOut_t arg, vectorOut_t darg);

        ErrorNormBased::ErrorNormBased(value_type alphaMin, value_type _a,
                                       value_type _b)
          : C (0.5 + alphaMin / 2), K ((1 - alphaMin) / 2), a (_a), b (_b)
        {}

        ErrorNormBased::ErrorNormBased(value_type alphaMin)
          : C (0.5 + alphaMin / 2), K ((1 - alphaMin) / 2)
        {
          static const value_type delta = 0.02;
          static const value_type r_half = 1e6;

          a = atanh ( (delta - 1 + C) / K ) / (1 - r_half);
          b = - r_half * a;
        }

        template bool ErrorNormBased::operator()
          (const HierarchicalIterative& solver, vectorOut_t arg, vectorOut_t darg);
      }

      HierarchicalIterative::HierarchicalIterative
      (const LiegroupSpacePtr_t& configSpace) :
        squaredErrorThreshold_ (0), inequalityThreshold_ (0),
        maxIterations_ (0), stacks_ (), configSpace_ (configSpace),
        dimension_ (0), reducedDimension_ (0), lastIsOptional_ (false),
        freeVariables_ (), saturate_ (), constraints_ (), lockedJoints_ (),
        sigma_ (0), dq_ (), dqSmall_ (), projector_ (), reducedJ_ (),
        saturation_ (configSpace->nv ()), reducedSaturation_ (),
        qSat_ (configSpace_->nq ()), tmpSat_ (), squaredNorm_ (0), datas_(),
        svd_ (), OM_ (configSpace->nv ()), OP_ (configSpace->nv ()),
        statistics_ ("HierarchicalIterative")
      {
        // Initialize freeVariables_ to all indices.
        freeVariables_.addRow (0, configSpace_->nv ());
      }

      HierarchicalIterative::HierarchicalIterative
      (const HierarchicalIterative& other) :
        squaredErrorThreshold_ (other.squaredErrorThreshold_),
        inequalityThreshold_ (other.inequalityThreshold_),
        maxIterations_ (other.maxIterations_), stacks_ (other.stacks_),
        configSpace_ (other.configSpace_), dimension_ (other.dimension_),
        reducedDimension_ (other.reducedDimension_),
        lastIsOptional_ (other.lastIsOptional_),
        freeVariables_ (other.freeVariables_),
        saturate_ (other.saturate_), constraints_ (other.constraints_),
        lockedJoints_ (other.lockedJoints_), sigma_(other.sigma_),
        dq_ (other.dq_), dqSmall_ (other.dqSmall_),
        projector_ (other.projector_),reducedJ_ (other.reducedJ_),
        saturation_ (other.saturation_),
        reducedSaturation_ (other.reducedSaturation_), qSat_ (other.qSat_),
        tmpSat_ (other.tmpSat_), squaredNorm_ (other.squaredNorm_),
        datas_ (other.datas_), svd_ (other.svd_), OM_ (other.OM_),
	OP_ (other.OP_), statistics_ (other.statistics_)
      {
      }

      bool HierarchicalIterative::contains
      (const ImplicitPtr_t& numericalConstraint) const
      {
        for (NumericalConstraints_t::const_iterator it = constraints_.begin ();
             it != constraints_.end (); ++it) {
          if (numericalConstraint == *it || *numericalConstraint == **it)
            return true;
        }
        return false;
      }

      void HierarchicalIterative::add (const DifferentiableFunctionPtr_t& f,
                                       const std::size_t& priority,
                                       const ComparisonTypes_t& comp)
      {
        add (Implicit::create (f, comp), priority);
      }

      void HierarchicalIterative::add (const ImplicitPtr_t& constraint,
                                       const std::size_t& priority)
      {
        const ComparisonTypes_t comp (constraint->comparisonType ());
        const std::size_t minSize = priority + 1;
        if (stacks_.size() < minSize) {
          stacks_.resize (minSize, ImplicitConstraintSet ());
        }
        stacks_ [priority].add (constraint);
        Data& d = datas_[priority];
        for (std::size_t i = 0; i < comp.size(); ++i) {
          switch (comp[i]) {
          case Superior:
          case Inferior:
            d.inequalityIndices.push_back (d.comparison.size());
            break;
          case Equality:
            d.equalityIndices.addRow(d.comparison.size(), 1);
            break;
          default:
            break;
          }
          d.comparison.push_back (comp[i]);
        }
        d.equalityIndices.updateRows<true, true, true>();
        update();

      }

      ArrayXb HierarchicalIterative::activeParameters () const
      {
        ArrayXb ap (ArrayXb::Constant(configSpace_->nq (), false));
        for (std::size_t i = 0; i < stacks_.size (); ++i) {
#ifndef NDEBUG
          dynamic_cast <const DifferentiableFunctionSet&>
            (stacks_[i].function ());
#endif
          const DifferentiableFunctionSet& dfs
            (dynamic_cast <const DifferentiableFunctionSet&>
             (stacks_[i].function ()));
          ap = ap || dfs.activeParameters();
        }
        return ap;
      }

      ArrayXb HierarchicalIterative::activeDerivativeParameters () const
      {
        ArrayXb ap (ArrayXb::Constant(configSpace_->nv (), false));
        for (std::size_t i = 0; i < stacks_.size (); ++i) {
#ifndef NDEBUG
          dynamic_cast <const DifferentiableFunctionSet&>
            (stacks_[i].function ());
#endif
          const DifferentiableFunctionSet& dfs
            (dynamic_cast <const DifferentiableFunctionSet&>
             (stacks_[i].function ()));
          ap = ap || dfs.activeDerivativeParameters();
        }
        return ap;
      }

      void HierarchicalIterative::update()
      {
        // Compute reduced size
        std::size_t reducedSize = freeVariables_.nbIndices();

        dimension_ = 0;
        reducedDimension_ = 0;
        for (std::size_t i = 0; i < stacks_.size (); ++i) {
          computeActiveRowsOfJ (i);

          const ImplicitConstraintSet& constraints (stacks_ [i]);
#ifndef NDEBUG
          dynamic_cast <const DifferentiableFunctionSet&>
            (constraints.function ());
#endif
          const DifferentiableFunctionSet& f
            (static_cast <const DifferentiableFunctionSet&>
             (constraints.function ()));
          dimension_ += f.outputSize();
          reducedDimension_ += datas_[i].activeRowsOfJ.nbRows();
          datas_[i].output = LiegroupElement (f.outputSpace ());
          datas_[i].rightHandSide = LiegroupElement (f.outputSpace ());
          datas_[i].rightHandSide.setNeutral ();

          assert(configSpace_->nv () == f.inputDerivativeSize());
          datas_[i].jacobian.resize(f.outputDerivativeSize(),
                                    f.inputDerivativeSize());
          datas_[i].jacobian.setZero();
          datas_[i].reducedJ.resize(datas_[i].activeRowsOfJ.nbRows(), reducedSize);

          datas_[i].svd = SVD_t (f.outputDerivativeSize(), reducedSize,
                                 Eigen::ComputeThinU | Eigen::ComputeThinV);
          datas_[i].svd.setThreshold (SVD_THRESHOLD);
          datas_[i].PK.resize (reducedSize, reducedSize);

          datas_[i].maxRank = 0;
        }

        dq_ = vector_t::Zero(configSpace_->nv ());
        dqSmall_.resize(reducedSize);
        projector_.resize(reducedSize, reducedSize);
        reducedJ_.resize(reducedDimension_, reducedSize);
        svd_ = SVD_t (reducedDimension_, reducedSize,
                      Eigen::ComputeThinU | Eigen::ComputeThinV);
      }

      void HierarchicalIterative::computeActiveRowsOfJ (std::size_t iStack)
      {
        Data& d = datas_[iStack];
        const ImplicitConstraintSet::Implicits_t constraints
          (stacks_ [iStack].constraints ());
        std::size_t row = 0;

        typedef Eigen::MatrixBlocks<false, false> BlockIndices;
        BlockIndices::segments_t rows;
        // Loop over functions of the stack
        for (std::size_t i = 0; i < constraints.size (); ++i) {
          ArrayXb adp = freeVariables_.rview
            (constraints [i]->function ().activeDerivativeParameters().
             matrix()).eval();
          if (adp.any()) // If at least one element of adp is true
            rows.push_back (BlockIndices::segment_t
                            (row, constraints [i]->function ().
                             outputDerivativeSize()));
          row += constraints [i]->function ().outputDerivativeSize();
        }
        d.activeRowsOfJ = Eigen::MatrixBlocks<false,false>
          (rows, freeVariables_.m_rows);
        d.activeRowsOfJ.updateRows<true, true, true>();
      }

      vector_t HierarchicalIterative::rightHandSideFromConfig
      (vectorIn_t config)
      {
        for (std::size_t i = 0; i < stacks_.size (); ++i) {
          ImplicitConstraintSet& ics = stacks_[i];
          Data& d = datas_[i];
          ics.function ().value (d.output, config);
          d.equalityIndices.lview(d.rightHandSide.vector ()) =
            d.equalityIndices.rview(d.output.vector ());
        }
        return rightHandSide();
      }

      bool HierarchicalIterative::rightHandSideFromConfig
      (const ImplicitPtr_t& constraint, vectorIn_t config)
      {
        const DifferentiableFunctionPtr_t& f (constraint->functionPtr ());
        for (std::size_t i = 0; i < stacks_.size (); ++i) {
          Data& d = datas_[i];
          const ImplicitConstraintSet& ics (stacks_[i]);
          assert (HPP_DYNAMIC_PTR_CAST (DifferentiableFunctionSet,
                                        ics.functionPtr ()));
          DifferentiableFunctionSetPtr_t dfs
            (HPP_STATIC_PTR_CAST (DifferentiableFunctionSet,
                                  ics.functionPtr ()));
          const DifferentiableFunctionSet::Functions_t& fs (dfs->functions ());
          size_type row = 0;
          for (std::size_t j = 0; j < fs.size(); ++j) {
            if (f == fs[j]) {
              LiegroupElement tmp (f->outputSpace ());
              f->value (tmp, config);
              d.output.vector ().segment(row, f->outputSize()) = tmp.vector ();
              for (size_type k = 0; k < f->outputSize(); ++k) {
                if (d.comparison[row + k] == Equality) {
                  d.rightHandSide.vector () [row + k] = d.output.vector ()[row + k];
                }
              }
              return true;
            }
            row += fs[j]->outputSize();
          }
        }
        return false;
      }

      bool HierarchicalIterative::rightHandSide
      (const ImplicitPtr_t& constraint, vectorIn_t rhs)
      {
        const DifferentiableFunctionPtr_t& f (constraint->functionPtr ());
        for (std::size_t i = 0; i < stacks_.size (); ++i) {
          Data& d = datas_[i];
          const ImplicitConstraintSet& ics (stacks_[i]);
          assert (HPP_DYNAMIC_PTR_CAST (DifferentiableFunctionSet,
                                        ics.functionPtr ()));
          DifferentiableFunctionSetPtr_t dfs
            (HPP_STATIC_PTR_CAST (DifferentiableFunctionSet,
                                  ics.functionPtr ()));
          const DifferentiableFunctionSet::Functions_t& fs (dfs->functions ());
          size_type row = 0;
          for (std::size_t j = 0; j < fs.size(); ++j) {
            if (f == fs[j]) {
              for (size_type k = 0; k < f->outputSize(); ++k) {
                if (d.comparison[row + k] == Equality) {
                  d.rightHandSide.vector () [row + k] = rhs [k];
                }
              }
              return true;
            }
            row += fs[j]->outputSize();
          }
        }
        return false;
      }

      void HierarchicalIterative::rightHandSide (vectorIn_t rhs)
      {
        size_type row = 0;
        for (std::size_t i = 0; i < stacks_.size (); ++i) {
          Data& d = datas_[i];
          d.equalityIndices.lview(d.rightHandSide.vector ())
            = rhs.segment(row, d.equalityIndices.m_nbRows);
          row += d.equalityIndices.m_nbRows;
        }
        assert (row == rhs.size());
      }

      vector_t HierarchicalIterative::rightHandSide () const
      {
        vector_t rhs(rightHandSideSize());
        size_type row = 0;
        for (std::size_t i = 0; i < stacks_.size (); ++i) {
          const Data& d = datas_[i];
          const size_type nRows = d.equalityIndices.m_nbRows;
          vector_t::SegmentReturnType seg = rhs.segment(row, nRows);
          seg = d.equalityIndices.rview(d.rightHandSide.vector ());
          row += nRows;
        }
        assert (row == rhs.size());
        return rhs;
      }

      size_type HierarchicalIterative::rightHandSideSize () const
      {
        size_type rhsSize = 0;
        for (std::size_t i = 0; i < stacks_.size (); ++i)
          rhsSize += datas_[i].equalityIndices.nbRows();
        return rhsSize;
      }

      template <bool ComputeJac>
      void HierarchicalIterative::computeValue (vectorIn_t config) const
      {
        for (std::size_t i = 0; i < stacks_.size (); ++i) {
          const ImplicitConstraintSet& constraints (stacks_ [i]);
          const DifferentiableFunction& f = constraints.function ();
          Data& d = datas_[i];

          f.value   (d.output, config);
          if (ComputeJac) f.jacobian(d.jacobian, config);
          d.error = d.output - d.rightHandSide;
          applyComparison<ComputeJac>(d.comparison, d.inequalityIndices,
                                      d.error, d.jacobian, inequalityThreshold_);

          // Copy columns that are not reduced
          if (ComputeJac) d.reducedJ = d.activeRowsOfJ.rview (d.jacobian);
        }
      }

      template void HierarchicalIterative::computeValue<false>(vectorIn_t config) const;
      template void HierarchicalIterative::computeValue<true >(vectorIn_t config) const;

      void HierarchicalIterative::computeSaturation (vectorIn_t config) const
      {
        bool applySaturate;
        applySaturate = saturate_ (config, qSat_, saturation_);
        if (!applySaturate) return;

        reducedSaturation_ = freeVariables_.rview (saturation_);
        assert (
                (    reducedSaturation_.array() == -1
                     || reducedSaturation_.array() ==  0
                     || reducedSaturation_.array() ==  1
                     ).all() );

        for (std::size_t i = 0; i < stacks_.size (); ++i) {
          Data& d = datas_[i];

          vector_t error = d.activeRowsOfJ.keepRows().rview(d.error);
          tmpSat_ = (reducedSaturation_.cast<value_type>().cwiseProduct
                     (d.reducedJ.transpose() * error).array() < 0);
          for (size_type j = 0; j < tmpSat_.size(); ++j)
            if (tmpSat_[j])
              d.reducedJ.col(j).setZero();
        }
      }

      void HierarchicalIterative::getValue (vectorOut_t v) const
      {
        size_type row = 0;
        for (std::size_t i = 0; i < datas_.size(); ++i) {
          const Data& d = datas_[i];
          v.segment(row, d.output.vector ().rows()) = d.output.vector ();
          row += d.output.vector ().rows();
        }
        assert (v.rows() == row);
      }

      void HierarchicalIterative::getReducedJacobian (matrixOut_t J) const
      {
        size_type row = 0;
        for (std::size_t i = 0; i < datas_.size(); ++i) {
          const Data& d = datas_[i];
          J.middleRows(row, d.reducedJ.rows()) = d.reducedJ;
          row += d.reducedJ.rows();
        }
        assert (J.rows() == row);
      }

      void HierarchicalIterative::computeError () const
      {
        const std::size_t end = (lastIsOptional_ ? stacks_.size() - 1 :
                                 stacks_.size());
        squaredNorm_ = 0;
        for (std::size_t i = 0; i < end; ++i) {
          const ImplicitConstraintSet::Implicits_t constraints
            (stacks_ [i].constraints ());
          const Data& d = datas_[i];
          size_type row = 0;
          for (std::size_t j = 0; j < constraints.size(); ++j) {
            size_type outputSize (constraints [j]->function ().outputSize ());
            squaredNorm_ = std::max
              (squaredNorm_,
               d.error.segment(row, outputSize).squaredNorm());
            row += outputSize;
          }
        }
      }

      void HierarchicalIterative::integrate
      (vectorIn_t from, vectorIn_t velocity, vectorOut_t result) const
      {
        typedef pinocchio::LiegroupElement Lge_t;
        typedef pinocchio::LiegroupConstElementRef LgeConstRef_t;
        LgeConstRef_t O (from, configSpace_);
        Lge_t M (O + velocity);
        saturate_ (M.vector (), result, saturation_);
      }

      void HierarchicalIterative::residualError (vectorOut_t error) const
      {
        size_type row = 0;
        for (std::size_t i = 0; i < datas_.size(); ++i) {
          const Data& d = datas_[i];
          error.segment(row, d.error.size()) = d.error;
          row += d.error.size();
        }
      }

      void HierarchicalIterative::computeDescentDirection () const
      {
        sigma_ = std::numeric_limits<value_type>::max();

        if (stacks_.empty()) {
          dq_.setZero();
          return;
        }
        vector_t err;
        if (stacks_.size() == 1) { // one level only
          Data& d = datas_[0];
          d.svd.compute (d.reducedJ);
          HPP_DEBUG_SVDCHECK (d.svd);
          // TODO Eigen::JacobiSVD does a dynamic allocation here.
          err = d.activeRowsOfJ.keepRows().rview(- d.error);
          dqSmall_ = d.svd.solve (err);
          d.maxRank = std::max(d.maxRank, d.svd.rank());
          if (d.maxRank > 0)
            sigma_ = std::min(sigma_, d.svd.singularValues()[d.maxRank - 1]);
        } else {
          projector_.setIdentity();
          for (std::size_t i = 0; i < stacks_.size (); ++i) {
            const DifferentiableFunction& f = stacks_[i].function ();
            Data& d = datas_[i];

            // TODO: handle case where this is the first element of the stack and it
            // has no functions
            if (f.outputSize () == 0) continue;
            /// projector is of size numberDof
            bool first = (i == 0);
            bool last = (i == stacks_.size() - 1);
            err = d.activeRowsOfJ.keepRows().rview(- d.error);
            if (first) {
              // dq should be zero and projector should be identity
              d.svd.compute (d.reducedJ);
              HPP_DEBUG_SVDCHECK (d.svd);
              // TODO Eigen::JacobiSVD does a dynamic allocation here.
              dqSmall_ = d.svd.solve (err);
            } else {
              // TODO check whether d.reducedJ * projector_ is not the null matrix.
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

            /// compute projector for next step.
            projectorOnSpan <SVD_t> (d.svd, d.PK);
            projector_ -= d.PK;
          }
        }
        expandDqSmall();
      }

      void HierarchicalIterative::expandDqSmall () const
      {
        Eigen::MatrixBlockView<vector_t, Eigen::Dynamic, 1, false, true>
          (dq_, freeVariables_.nbIndices(), freeVariables_.indices()) =
          dqSmall_;
      }

      std::ostream& HierarchicalIterative::print (std::ostream& os) const
      {
        os << "HierarchicalIterative, " << stacks_.size() << " level." << iendl
           << "dimension " << dimension() << iendl
           << "reduced dimension " << reducedDimension() << iendl
           << "reduction: " << freeVariables_ << incindent;
        const std::size_t end = (lastIsOptional_ ? stacks_.size() - 1 :
                                 stacks_.size());
        for (std::size_t i = 0; i < stacks_.size(); ++i) {
          const ImplicitConstraintSet::Implicits_t constraints
            (stacks_ [i].constraints ());
          const Data& d = datas_[i];
          os << iendl << "Level " << i;
          if (lastIsOptional_ && i == end) os << " (optional)";
          os << ": Stack of " << constraints.size () << " functions" << incindent;
          size_type row = 0;
          for (std::size_t j = 0; j < constraints.size (); ++j) {
            const DifferentiableFunctionPtr_t& f
              (constraints [j]->functionPtr ());
            os << iendl << j << ": ["
               << row << ", " << f->outputSize() << "],"
               << *f
               << iendl << "Rhs: " << condensed(d.rightHandSide.vector().segment
                                                (row, f->outputSize()));
            row += f->outputSize();
          }
          os << decendl;
          os << "Equality idx: " << d.equalityIndices;
          os << iendl << "Active rows: " << d.activeRowsOfJ;
        }
        return os << decindent;
      }

      template HierarchicalIterative::Status HierarchicalIterative::solve
      (vectorOut_t arg, lineSearch::Constant       lineSearch) const;
      template HierarchicalIterative::Status HierarchicalIterative::solve
      (vectorOut_t arg, lineSearch::Backtracking   lineSearch) const;
      template HierarchicalIterative::Status HierarchicalIterative::solve
      (vectorOut_t arg, lineSearch::FixedSequence  lineSearch) const;
      template HierarchicalIterative::Status HierarchicalIterative::solve
      (vectorOut_t arg, lineSearch::ErrorNormBased lineSearch) const;
    } // namespace solver
  } // namespace constraints
} // namespace hpp
