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

#include <hpp/constraints/explicit-solver.hh>

#include <queue>

#include <hpp/util/indent.hh>

#include <hpp/pinocchio/util.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/liegroup.hh>

#include <hpp/constraints/matrix-view.hh>


namespace hpp {
  namespace constraints {
    typedef Eigen::MatrixBlocksRef<false, false> MatrixBlocksRef;

    namespace {
      void append (const Eigen::RowBlockIndices& rbi, std::queue<size_type>& q) {
        for (std::size_t i = 0; i < rbi.indices().size(); ++i)
          for (size_type j = 0; j < rbi.indices()[i].second; ++j)
            q.push(rbi.indices()[i].first + j);
      }
    }

    Eigen::ColBlockIndices ExplicitSolver::activeParameters () const
    {
      return inArgs_.transpose();
    }

    const Eigen::ColBlockIndices& ExplicitSolver::activeDerivativeParameters () const
    {
      return inDers_;
    }

    bool ExplicitSolver::solve (vectorOut_t arg) const
    {
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        computeFunction(computationOrder_[i], arg);
      }
      return true;
    }

    bool ExplicitSolver::isSatisfied (vectorIn_t arg, vectorOut_t error) const
    {
      value_type squaredNorm = 0;

      size_type row = 0;
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        const Function& f = functions_[i];
        // Compute this function
        f.qin = f.inArg.rview(arg);
        f.f->value(f.value, f.qin);
        f.value += f.rightHandSide;
        const size_type& nbRows = f.outDer.nbRows();
        f.qout = f.outArg.rview(arg);
        if (f.g) f.g->value(f.expected, f.qout);
        else     f.expected.vector() = f.qout;
        error.segment (row, nbRows) = f.expected - f.value;
        squaredNorm = std::max(squaredNorm,
            error.segment (row, nbRows).squaredNorm ());
        row += nbRows;
      }
      assert (row == error.size());
      hppDout (info, "Max squared error norm is " << squaredNorm);
      return squaredNorm < squaredErrorThreshold_;
    }

    bool ExplicitSolver::isSatisfied (vectorIn_t arg) const
    {
      diffSmall_.resize(outDers_.nbIndices());
      return isSatisfied (arg, diffSmall_);
    }

    ExplicitSolver::Function::Function (DifferentiableFunctionPtr_t _f,
        RowBlockIndices ia, RowBlockIndices oa,
        ColBlockIndices id, RowBlockIndices od,
        const ComparisonTypes_t& comp) :
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

    void ExplicitSolver::Function::setG (const DifferentiableFunctionPtr_t& _g,
                                         const DifferentiableFunctionPtr_t& _ginv)
    {
      assert ( (!_g && !_ginv) // No function g and ginv
          || ( _g && _ginv
            && * f->outputSpace() == *_g->outputSpace()
            && *_g->outputSpace() == *_ginv->outputSpace())
          // Functions specified and the output spaces are all the same.
          );
      g = _g;
      ginv = _ginv;
      size_type n = (g ? g->outputSpace()->nv() : 0);
      jGinv.resize (n,n);
    }

    size_type ExplicitSolver::add (const DifferentiableFunctionPtr_t& f,
        const RowBlockIndices& inArg,
        const RowBlockIndices& outArg,
        const ColBlockIndices& inDer,
        const RowBlockIndices& outDer)
    {
      return add (f, inArg, outArg, inDer, outDer,
          ComparisonTypes_t(f->outputDerivativeSize(), EqualToZero));
    }

    size_type ExplicitSolver::add (const DifferentiableFunctionPtr_t& f,
        const RowBlockIndices& inArg,
        const RowBlockIndices& outArg,
        const ColBlockIndices& inDer,
        const RowBlockIndices& outDer,
        const ComparisonTypes_t& comp)
    {
      assert (outArg.indices().size() == 1 && "Only contiguous function output is supported.");
      assert (outDer.indices().size() == 1 && "Only contiguous function output is supported.");
      const RowBlockIndices::segment_t& outIdx = outArg.indices()[0];
      const RowBlockIndices::segment_t& outDerIdx = outDer.indices()[0];

      // TODO This sanity check should not be necessary
      // It is done in the while loop below
      // Sanity check: is it explicit ?
      for (std::size_t i = 0; i < inArg.indices().size(); ++i)
        if (BlockIndex::overlap(inArg.indices()[i], outIdx))
          return -1;
      // Sanity check: Comparison type must be either EqualToZero or Equality
      assert (comp.size() == (std::size_t)f->outputDerivativeSize());
      for (std::size_t i = 0; i < comp.size(); ++i)
        if (comp[i] != EqualToZero && comp[i] != Equality) return -1;
      // Check that no other function already computes its outputs.
      if ((outArg.rview(argFunction_).eval().array() >= 0).any())
        return -1;
      // Check that it does not insert a loop.
      std::queue<size_type> idxArg;
      append(inArg, idxArg);
      while (!idxArg.empty()) {
        // iArg must be computed before f
        size_type iArg = idxArg.back();
        idxArg.pop();
        // iArg is an output of f -> cannot be computed before f
        if (iArg >= outIdx.first && iArg < outIdx.first + outIdx.second) return -1;
        // iArg is not computed by any function
        if (argFunction_[iArg] < 0) continue;
        const Function& func = functions_[argFunction_[iArg]];
        append(func.inArg, idxArg);
      }

      // Add the function
      int idx = int(functions_.size());
      outArg.lview(argFunction_).setConstant(idx);
      outDer.lview(derFunction_).setConstant(idx);
      functions_.push_back (Function(f, inArg, outArg, inDer, outDer, comp));

      // Update the free dofs
      outArgs_.addRow(outIdx.first, outIdx.second);
      outArgs_.updateIndices<true, true, true>();
      freeArgs_ = RowBlockIndices
        (BlockIndex::difference (BlockIndex::segment_t(0, argSize_),
                                 outArgs_.indices()));

      BlockIndex::add (inArgs_.m_rows, inArg.rows());
      inArgs_ = RowBlockIndices
        (BlockIndex::difference (inArgs_.rows(), outArgs_.rows()));
      // should be sorted already
      inArgs_.updateIndices<false, true, true>();

      outDers_.addRow(outDerIdx.first, outDerIdx.second);
      outDers_.updateIndices<true, true, true>();
      freeDers_ = ColBlockIndices
        (BlockIndex::difference(BlockIndex::segment_t(0, derSize_),
                                outDers_.indices()));

      BlockIndex::add (inDers_.m_cols, inDer.cols());
      inDers_ = ColBlockIndices
        (BlockIndex::difference (inDers_.cols(), outDers_.rows()));
      // should be sorted already
      inDers_.updateIndices<false, true, true>();

      /// Computation order
      std::size_t order = 0;
      computationOrder_.resize(functions_.size());
      inOutDependencies_ = Eigen::MatrixXi::Zero(functions_.size(), derSize_);
      Computed_t computed(functions_.size(), false);
      for(std::size_t i = 0; i < functions_.size(); ++i)
        computeOrder(i, order, computed);
      assert(order == functions_.size());
      return functions_.size() - 1;
    }

    bool ExplicitSolver::setG (const DifferentiableFunctionPtr_t& df,
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

    bool ExplicitSolver::replace (
        const DifferentiableFunctionPtr_t& oldf,
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

    void ExplicitSolver::computeFunction(const std::size_t& iF, vectorOut_t arg) const
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

    void ExplicitSolver::jacobian(matrixOut_t jacobian, vectorIn_t arg) const
    {
      // TODO this could be done only on the complement of inDers_
      jacobian.setZero();
      MatrixBlocksRef (freeDers_, freeDers_)
        .lview (jacobian).setIdentity();
      // Compute the function jacobians
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        const Function& f = functions_[i];
        f.qin = f.inArg.rview(arg);
        if (f.ginv) f.f->value(f.value, f.qin);
        f.f->jacobian(f.jacobian, f.qin);
        if (f.equalityIndices.nbIndices() > 0)
          f.f->outputSpace ()->Jintegrate (f.rightHandSide, f.jacobian);
        if (f.ginv) {
          f.value += f.rightHandSide;
          f.ginv->jacobian(f.jGinv, f.value.vector());
          f.jacobian.applyOnTheLeft(f.jGinv);
        }
      }
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        computeJacobian(computationOrder_[i], jacobian);
      }
    }

    void ExplicitSolver::computeJacobian(const std::size_t& iF, matrixOut_t J) const
    {
      const Function& f = functions_[iF];
      matrix_t Jg (MatrixBlocksRef (f.inDer, inDers_).rview(J));
      MatrixBlocksRef (f.outDer, inDers_).lview (J) = f.jacobian * Jg;
    }

    void ExplicitSolver::computeOrder(const std::size_t& iF, std::size_t& iOrder, Computed_t& computed)
    {
      if (computed[iF]) return;
      const Function& f = functions_[iF];
      for (std::size_t i = 0; i < f.inDer.indices().size(); ++i) {
        const BlockIndex::segment_t& segment = f.inDer.indices()[i];
        for (size_type j = 0; j < segment.second; ++j) {
          if (derFunction_[segment.first + j] < 0) {
            inOutDependencies_(iF, segment.first + j) += 1;
          } else {
            assert((std::size_t)derFunction_[segment.first + j] < functions_.size());
            computeOrder(derFunction_[segment.first + j], iOrder, computed);
            inOutDependencies_.row(iF) += inOutDependencies_.row(derFunction_[segment.first + j]);
          }
        }
      }
      computationOrder_[iOrder] = iF;
      ++iOrder;
      computed[iF] = true;
    }

    vector_t ExplicitSolver::rightHandSideFromInput (vectorIn_t arg)
    {
      for (std::size_t i = 0; i < functions_.size (); ++i) {
        Function& f = functions_[i];
        f.qin = f.inArg.rview(arg);
        f.f->value(f.value, f.qin);
        f.qout = f.outArg.rview(arg);
        if (f.g) f.g->value(f.expected, f.qout);
        else     f.expected.vector() = f.qout;
        vector_t rhs = f.expected - f.value;
        f.equalityIndices.lview(f.rightHandSide) = f.equalityIndices.rview(rhs);
      }
      return rightHandSide();
    }

    bool ExplicitSolver::rightHandSideFromInput (const DifferentiableFunctionPtr_t& df, vectorIn_t arg)
    {
      for (std::size_t i = 0; i < functions_.size (); ++i) {
        Function& f = functions_[i];
        if (f.f == df) {
          rightHandSideFromInput (i, arg);
          return true;
        }
      }
      return false;
    }

    void ExplicitSolver::rightHandSideFromInput (const size_type& fidx, vectorIn_t arg)
    {
      assert (fidx < functions_.size());
      Function& f = functions_[fidx];

      // Computes f(q1) and g(q2)
      f.qin = f.inArg.rview(arg);
      f.f->value(f.value, f.qin);
      f.qout = f.outArg.rview(arg);
      if (f.g) f.g->value(f.expected, f.qout);
      else     f.expected.vector() = f.qout;

      // Set rhs = g(q2) - f(q1)
      vector_t rhs = f.expected - f.value;
      f.equalityIndices.lview(f.rightHandSide) = f.equalityIndices.rview(rhs);
    }

    bool ExplicitSolver::rightHandSide (const DifferentiableFunctionPtr_t& df, vectorIn_t rhs)
    {
      for (std::size_t i = 0; i < functions_.size (); ++i) {
        Function& f = functions_[i];
        if (f.f == df) {
          rightHandSide (i, rhs);
          return true;
        }
      }
      return false;
    }

    void ExplicitSolver::rightHandSide (const size_type& i, vectorIn_t rhs)
    {
      assert (i < functions_.size());
      Function& f = functions_[i];
      f.equalityIndices.lview(f.rightHandSide) = f.equalityIndices.rview(rhs);
    }

    void ExplicitSolver::rightHandSide (vectorIn_t rhs)
    {
      size_type row = 0;
      for (std::size_t i = 0; i < functions_.size (); ++i) {
        Function& f = functions_[i];

        f.equalityIndices.lview(f.rightHandSide)
          = rhs.segment(row, f.equalityIndices.nbRows());
        row += f.equalityIndices.nbRows();
      }
      assert (row == rhs.size());
    }

    vector_t ExplicitSolver::rightHandSide () const
    {
      vector_t rhs(rightHandSideSize());
      size_type row = 0;
      for (std::size_t i = 0; i < functions_.size (); ++i) {
        const Function& f = functions_[i];
        const size_type nRows = f.equalityIndices.nbRows();
        vector_t::SegmentReturnType seg = rhs.segment(row, nRows);
        seg = f.equalityIndices.rview(f.rightHandSide);
        row += nRows;
      }
      assert (row == rhs.size());
      return rhs;
    }

    size_type ExplicitSolver::rightHandSideSize () const
    {
      size_type rhsSize = 0;
      for (std::size_t i = 0; i < functions_.size (); ++i)
        rhsSize += functions_[i].equalityIndices.nbRows();
      return rhsSize;
    }

    std::ostream& ExplicitSolver::print (std::ostream& os) const
    {
      os << "ExplicitSolver, " << functions_.size() << " functions." << incendl
        << "Free args: " << freeArgs_ << iendl
        << "Params: " << inArgs_ << " -> " << outArgs_ << iendl
        << "Dofs: "   << inDers_ << " -> " << outDers_ << iendl
        << "Functions" << incindent;
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        const Function& f = functions_[computationOrder_[i]];
        os << iendl << i << ": " << f.inArg << " -> " << f.outArg
          << incendl << *f.f
          << decendl << "Rhs: " << condensed(f.rightHandSide)
          << iendl   << "Equality: " << f.equalityIndices;
      }
      return os << decindent << decindent;
    }

    Eigen::MatrixXi ExplicitSolver::inOutDofDependencies () const
    {
      Eigen::MatrixXi iod (derSize(), inDers_.nbCols());
      if (inDers_.nbCols() == 0) return iod;
      Eigen::RowVectorXi tmp (inDers_.nbCols());
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        const Function& f = functions_[i];
        tmp = inDers_.rview (inOutDependencies_.row(i));
        f.outDer.lview(iod).rowwise() = tmp;
      }
      return outDers_.rview(iod);
    }
  } // namespace constraints
} // namespace hpp
