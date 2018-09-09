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

#include <hpp/constraints/explicit-constraint-set.hh>

#include <queue>

#include <hpp/util/indent.hh>

#include <hpp/pinocchio/util.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/liegroup.hh>

#include <hpp/constraints/matrix-view.hh>
#include <hpp/constraints/explicit.hh>


namespace hpp {
  namespace constraints {
    typedef Eigen::MatrixBlocksRef<false, false> MatrixBlocksRef;

    namespace {
      /// Append all indices of a set of segments to a queue of indices
      void append (const Eigen::RowBlockIndices& rbi, std::queue<size_type>& q)
      {
        append (rbi.indices (), q);
    }
      /// Append all indices of a set of segments to a queue of indices
      void append (const segments_t& segments, std::queue<size_type>& q)
      {
        for (std::size_t i = 0; i < segments.size(); ++i)
          for (size_type j = 0; j < segments[i].second; ++j)
            q.push (segments [i].first + j);
      }
    }

    Eigen::ColBlockIndices ExplicitConstraintSet::activeParameters () const
    {
      return inArgs_.transpose();
    }

    const Eigen::ColBlockIndices&
    ExplicitConstraintSet::activeDerivativeParameters () const
    {
      return inDers_;
    }

    bool ExplicitConstraintSet::solve (vectorOut_t arg) const
    {
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        computeFunction(computationOrder_[i], arg);
      }
      return true;
    }

    bool ExplicitConstraintSet::isSatisfied (vectorIn_t arg, vectorOut_t error)
      const
    {
      value_type squaredNorm = 0;

      size_type row = 0;
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        const Function& f = functions_[i];
        // Compute this function
        f.qin = f.inArg.rview(arg); // q_{in}
        f.f->value(f.value, f.qin); // f (q_{in})
        f.value += f.rightHandSide; // f (q_{in}) + rhs
        const size_type& nbRows = f.outDer.nbRows();
        f.qout = f.outArg.rview(arg); // q_{out}
        if (f.g) f.g->value(f.expected, f.qout); // g (q_{out})
        else     f.expected.vector() = f.qout;
        // g (q_{out}) - (f (q_{in}) + rhs)
        error.segment (row, nbRows) = f.expected - f.value;
        squaredNorm = std::max(squaredNorm,
            error.segment (row, nbRows).squaredNorm ());
        row += nbRows;
      }
      assert (row == error.size());
      hppDout (info, "Max squared error norm is " << squaredNorm);
      return squaredNorm < squaredErrorThreshold_;
    }

    bool ExplicitConstraintSet::isSatisfied (vectorIn_t arg) const
    {
      diffSmall_.resize(outDers_.nbIndices());
      return isSatisfied (arg, diffSmall_);
    }

    ExplicitConstraintSet::Function::Function
    (DifferentiableFunctionPtr_t _f, RowBlockIndices ia, RowBlockIndices oa,
     ColBlockIndices id, RowBlockIndices od, const ComparisonTypes_t& comp) :
      f (_f), inArg (ia), outArg (oa), inDer (id), outDer (od),
      comparison (comp),
      rightHandSide (vector_t::Zero(f->outputSpace()->nv())),
      value (f->outputSpace ()), expected (f->outputSpace ())
    {
      jacobian.resize(_f->outputDerivativeSize(), _f->inputDerivativeSize());
      for (std::size_t i = 0; i < comp.size(); ++i) {
        switch (comp[i]) {
          case Equality:
            equalityIndices.addRow(i, 1);
            break;
          case Superior:
          case Inferior:
          default:
            break;
        }
      }
      equalityIndices.updateRows<true, true, true>();
    }

    void ExplicitConstraintSet::Function::setG
    (const DifferentiableFunctionPtr_t& _g,
     const DifferentiableFunctionPtr_t& _ginv)
    size_type ExplicitConstraintSet::add (const ExplicitPtr_t& constraint)
    {
      assert (constraint->outputConf ().size() == 1 &&
              "Only contiguous function output is supported.");
      assert (constraint->outputVelocity ().size() == 1 &&
              "Only contiguous function output is supported.");
      typedef Eigen::RowBlockIndices RowBlockIndices;
      typedef Eigen::ColBlockIndices ColBlockIndices;
      const RowBlockIndices::segment_t& outIdx (constraint->outputConf () [0]);
      const RowBlockIndices::segment_t& outDerIdx
        (constraint->outputVelocity () [0]);

      // TODO This sanity check should not be necessary
      // It is done in the while loop below
      // Sanity check: is it explicit ?
      for (std::size_t i = 0; i < constraint->inputConf ().size(); ++i)
        if (BlockIndex::overlap(constraint->inputConf ()[i], outIdx)) {
          size_type f0 (constraint->inputConf ()[i].first);
          size_type s0 (f0 + constraint->inputConf ()[i].second - 1);
          size_type f1 (outIdx.first);
          size_type s1 (f1 + outIdx.second - 1);
          std::ostringstream oss;
          oss << "Explicit constraint \"" << constraint->function ().name ()
              << "\" is malformed.";
          oss << " input [" << f0 << "," << s0
              << "] and output [" << f1 << "," << s1 << "] segments overlap.";
            throw std::logic_error (oss.str ().c_str ());
        }
      // Sanity check: Comparison type must be either EqualToZero or Equality
      const ComparisonTypes_t& comp (constraint->comparisonType ());
      assert (comp.size() ==
              (std::size_t) constraint->function ().outputDerivativeSize());
      for (std::size_t i = 0; i < comp.size(); ++i)
        if (comp[i] != EqualToZero && comp[i] != Equality) return -1;
      // Check that no other function already computes its outputs.
      RowBlockIndices outArg (constraint->outputConf ());
      if ((outArg.rview(argFunction_).eval().array() >= 0).any())
        return -1;
      // Check that it does not insert a loop.
      std::queue<size_type> idxArg;
      // Put all indices of inArg into idxArg
      RowBlockIndices inArg (constraint->inputConf ());
      append(inArg, idxArg);
      while (!idxArg.empty()) {
        // iArg must be computed before f
        size_type iArg = idxArg.back();
        idxArg.pop();
        // iArg is an output of the new function -> new function should be
        // computed earlier -> incompatible.
        if (iArg >= outIdx.first && iArg < outIdx.first + outIdx.second)
          return -1;
        // iArg is not computed by any function
        if (argFunction_ [iArg] < 0) continue;
        const Data& d = data_ [argFunction_ [iArg]];
        append(d.constraint->inputConf (), idxArg);
      }

      // Add the function
      int idx = int(data_.size());
      outArg.lview(argFunction_).setConstant(idx);
      RowBlockIndices (constraint->outputVelocity ()).lview(derFunction_).
        setConstant(idx);
      data_.push_back (Data (constraint));

      // Update the free dofs
      outArgs_.addRow(outIdx.first, outIdx.second);
      outArgs_.updateIndices<true, true, true>();
      notOutArgs_ = RowBlockIndices
        (BlockIndex::difference (BlockIndex::segment_t(0, nq_),
                                 outArgs_.indices()));

      BlockIndex::add (inArgs_.m_rows, inArg.rows());
      inArgs_ = RowBlockIndices
        (BlockIndex::difference (inArgs_.rows(), outArgs_.rows()));
      // should be sorted already
      inArgs_.updateIndices<false, true, true>();

      outDers_.addRow(outDerIdx.first, outDerIdx.second);
      outDers_.updateIndices<true, true, true>();
      notOutDers_ = ColBlockIndices
        (BlockIndex::difference(BlockIndex::segment_t(0, nv_),
                                outDers_.indices()));

      BlockIndex::add (inDers_.m_cols,
                       ColBlockIndices (constraint->inputVelocity ()).cols());
      inDers_ = ColBlockIndices
        (BlockIndex::difference (inDers_.cols(), outDers_.rows()));
      // should be sorted already
      inDers_.updateIndices<false, true, true>();

      /// Computation order
      std::size_t order = 0;
      computationOrder_.resize(functions_.size());
      inOutDependencies_ = Eigen::MatrixXi::Zero(functions_.size(), nv_);
      Computed_t computed(functions_.size(), false);
      for(std::size_t i = 0; i < functions_.size(); ++i)
        computeOrder(i, order, computed);
      assert(order == functions_.size());
      return functions_.size() - 1;
    }

    bool ExplicitConstraintSet::setG (const DifferentiableFunctionPtr_t& df,
                                      const DifferentiableFunctionPtr_t& g,
                                      const DifferentiableFunctionPtr_t& ginv)
    {
      for (std::size_t i = 0; i < functions_.size (); ++i) {
        Function& f = functions_[i];
        if (f.f == df) {
          f.setG (g, ginv);
          return true;
        }
      }
      return false;
    }

    bool ExplicitConstraintSet::replace
    (const DifferentiableFunctionPtr_t& oldf,
     const DifferentiableFunctionPtr_t& newf)
    {
      assert(oldf->inputSize() == newf->inputSize()
          && oldf->inputDerivativeSize() == newf->inputDerivativeSize()
          && oldf->outputSize() == newf->outputSize()
          && oldf->outputDerivativeSize() == newf->outputDerivativeSize());
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        if (functions_[i].f == oldf) {
          functions_[i].f = newf;
          return true;
        }
      }
      return false;
    }

    void ExplicitConstraintSet::computeFunction
    (const std::size_t& iF, vectorOut_t arg) const
    {
      const Function& f = functions_[iF];
      // Compute this function
      f.qin = f.inArg.rview(arg);
      f.f->value(f.value, f.qin);
      f.value += f.rightHandSide;
      if (f.ginv) f.ginv->value (f.expected, f.value.vector());
      else        f.expected.vector() = f.value.vector();
      f.outArg.lview(arg) = f.expected.vector();
    }

    void ExplicitConstraintSet::jacobian
    (matrixOut_t jacobian, vectorIn_t arg) const
    {
      // TODO this could be done only on the complement of inDers_
      jacobian.setZero();
      MatrixBlocksRef (notOutDers_, notOutDers_)
        .lview (jacobian).setIdentity();
      // Compute the function jacobians
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        const Function& E = functions_[i];
        E.qin = E.inArg.rview(arg);
        if (E.ginv) E.f->value(E.value, E.qin);
        E.f->jacobian(E.jacobian, E.qin);
        if (E.equalityIndices.nbIndices() > 0)
          E.f->outputSpace ()->dIntegrate_dq (E.value, E.rightHandSide, E.jacobian);
        if (E.ginv) {
          E.value += E.rightHandSide;
          E.ginv->jacobian(E.jGinv, E.value.vector());
          E.jacobian.applyOnTheLeft(E.jGinv);
        }
      }
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        computeJacobian(computationOrder_[i], jacobian);
      }
    }

    void ExplicitConstraintSet::computeJacobian
    (const std::size_t& iE, matrixOut_t J) const
    {
      const Function& E = functions_[iE];
      matrix_t Jin (MatrixBlocksRef (E.inDer, inDers_).rview(J));
      // Jout = E.jacobian * Jin
      MatrixBlocksRef (E.outDer, inDers_).lview (J) = E.jacobian * Jin;
    }

    void ExplicitConstraintSet::computeOrder
    (const std::size_t& iE, std::size_t& iOrder, Computed_t& computed)
    {
      if (computed[iE]) return;
      const Function& E = functions_[iE];
      for (std::size_t i = 0; i < E.inDer.indices().size(); ++i) {
        const BlockIndex::segment_t& segment = E.inDer.indices()[i];
        for (size_type j = 0; j < segment.second; ++j) {
          if (derFunction_[segment.first + j] < 0) {
            inOutDependencies_(iE, segment.first + j) += 1;
          } else {
            assert((std::size_t)derFunction_[segment.first + j] < functions_.size());
            computeOrder(derFunction_[segment.first + j], iOrder, computed);
            inOutDependencies_.row(iE) += inOutDependencies_.row(derFunction_[segment.first + j]);
          }
        }
      }
      computationOrder_[iOrder] = iE;
      ++iOrder;
      computed[iE] = true;
    }

    vector_t ExplicitConstraintSet::rightHandSideFromInput (vectorIn_t arg)
    {
      for (std::size_t i = 0; i < functions_.size (); ++i) {
        Function& E = functions_[i];
        E.qin = E.inArg.rview(arg);
        E.f->value(E.value, E.qin);
        E.qout = E.outArg.rview(arg);
        if (E.g) E.g->value(E.expected, E.qout);
        else     E.expected.vector() = E.qout;
        vector_t rhs = E.expected - E.value;
        E.equalityIndices.lview(E.rightHandSide) = E.equalityIndices.rview(rhs);
      }
      return rightHandSide();
    }

    bool ExplicitConstraintSet::rightHandSideFromInput
    (const DifferentiableFunctionPtr_t& df, vectorIn_t arg)
    {
      for (std::size_t i = 0; i < functions_.size (); ++i) {
        Function& E = functions_[i];
        if (E.f == df) {
          rightHandSideFromInput (i, arg);
          return true;
        }
      }
      return false;
    }

    void ExplicitConstraintSet::rightHandSideFromInput
    (const size_type& fidx, vectorIn_t arg)
    {
      assert (fidx < (size_type) functions_.size());
      Function& E = functions_[fidx];

      // Computes f(q1) and g(q2)
      E.qin = E.inArg.rview(arg);
      E.f->value(E.value, E.qin);
      E.qout = E.outArg.rview(arg);
      if (E.g) E.g->value(E.expected, E.qout);
      else     E.expected.vector() = E.qout;

      // Set rhs = g(q2) - f(q1)
      vector_t rhs = E.expected - E.value;
      E.equalityIndices.lview(E.rightHandSide) = E.equalityIndices.rview(rhs);
    }

    bool ExplicitConstraintSet::rightHandSide
    (const DifferentiableFunctionPtr_t& df, vectorIn_t rhs)
    {
      for (std::size_t i = 0; i < functions_.size (); ++i) {
        Function& E = functions_[i];
        if (E.f == df) {
          rightHandSide (i, rhs);
          return true;
        }
      }
      return false;
    }

    void ExplicitConstraintSet::rightHandSide
    (const size_type& i, vectorIn_t rhs)
    {
      assert (i < (size_type) functions_.size());
      Function& E = functions_[i];
      E.equalityIndices.lview(E.rightHandSide) = E.equalityIndices.rview(rhs);
    }

    void ExplicitConstraintSet::rightHandSide (vectorIn_t rhs)
    {
      size_type row = 0;
      for (std::size_t i = 0; i < functions_.size (); ++i) {
        Function& E = functions_[i];

        E.equalityIndices.lview(E.rightHandSide)
          = rhs.segment(row, E.equalityIndices.nbRows());
        row += E.equalityIndices.nbRows();
      }
      assert (row == rhs.size());
    }

    vector_t ExplicitConstraintSet::rightHandSide () const
    {
      vector_t rhs(rightHandSideSize());
      size_type row = 0;
      for (std::size_t i = 0; i < functions_.size (); ++i) {
        const Function& E = functions_[i];
        const size_type nRows = E.equalityIndices.nbRows();
        vector_t::SegmentReturnType seg = rhs.segment(row, nRows);
        seg = E.equalityIndices.rview(E.rightHandSide);
        row += nRows;
      }
      assert (row == rhs.size());
      return rhs;
    }

    size_type ExplicitConstraintSet::rightHandSideSize () const
    {
      size_type rhsSize = 0;
      for (std::size_t i = 0; i < functions_.size (); ++i)
        rhsSize += functions_[i].equalityIndices.nbRows();
      return rhsSize;
    }

    std::ostream& ExplicitConstraintSet::print (std::ostream& os) const
    {
      os << "ExplicitConstraintSet, " << functions_.size()
         << " functions." << incendl
        << "Other args: " << notOutArgs_ << iendl
        << "Params: " << inArgs_ << " -> " << outArgs_ << iendl
        << "Dofs: "   << inDers_ << " -> " << outDers_ << iendl
        << "Functions" << incindent;
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        const Function& E = functions_[computationOrder_[i]];
        os << iendl << i << ": " << E.inArg << " -> " << E.outArg
          << incendl << *E.f
          << decendl << "Rhs: " << condensed(E.rightHandSide)
          << iendl   << "Equality: " << E.equalityIndices;
      }
      return os << decindent << decindent;
    }

    Eigen::MatrixXi ExplicitConstraintSet::inOutDofDependencies () const
    {
      Eigen::MatrixXi iod (nv (), inDers_.nbCols());
      if (inDers_.nbCols() == 0) return iod;
      Eigen::RowVectorXi tmp (inDers_.nbCols());
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        const Function& E = functions_[i];
        tmp = inDers_.rview (inOutDependencies_.row(i));
        E.outDer.lview(iod).rowwise() = tmp;
      }
      return outDers_.rview(iod);
    }
  } // namespace constraints
} // namespace hpp
