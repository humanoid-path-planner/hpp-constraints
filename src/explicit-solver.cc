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

namespace hpp {
  namespace constraints {
    typedef Eigen::BlockIndex<size_type> BlockIndex;
    typedef Eigen::MatrixBlockView<matrixOut_t, Eigen::Dynamic, Eigen::Dynamic, false, false> MatrixOutView_t;

    namespace {
      void append (const Eigen::RowBlockIndexes& rbi, std::queue<size_type>& q) {
        for (std::size_t i = 0; i < rbi.indexes().size(); ++i)
          for (size_type j = 0; j < rbi.indexes()[i].second; ++j)
            q.push(rbi.indexes()[i].first + j);
      }
    }

    bool ExplicitSolver::solve (vectorOut_t arg) const
    {
      Computed_t computed(functions_.size(), false);
      for(std::size_t i = 0; i < computed.size(); ++i) {
        computeFunction(i, arg, computed);
      }
      return true;
    }

    bool ExplicitSolver::add (const DifferentiableFunctionPtr_t& f,
        const RowBlockIndexes& inArg,
        const RowBlockIndexes& outArg,
        const ColBlockIndexes& inDer,
        const RowBlockIndexes& outDer)
    {
      assert (outArg.indexes().size() == 1 && "Only contiguous function output is supported.");
      assert (outDer.indexes().size() == 1 && "Only contiguous function output is supported.");
      const RowBlockIndexes::BlockIndexType& outIdx = outArg.indexes()[0];
      const RowBlockIndexes::BlockIndexType& outDerIdx = outDer.indexes()[0];

      // Sanity check: is it explicit ?
      for (std::size_t i = 0; i < inArg.indexes().size(); ++i)
        if (BlockIndex::overlap(inArg.indexes()[i], outIdx))
          return false;
      // Check that no other function already computes its outputs.
      if ((outArg.view(argFunction_).eval().array() >= 0).any())
        return false;
      // Check that it does not insert a loop.
      std::queue<size_type> idxArg;
      append(inArg, idxArg);
      while (!idxArg.empty()) {
        size_type iArg = idxArg.back();
        idxArg.pop();
        if (argFunction_[iArg] < 0) continue;
        if (iArg >= outIdx.first && iArg < outIdx.first + outIdx.second) return false;
        const Function& func = functions_[argFunction_[iArg]];
        append(func.inArg, idxArg);
      }

      // Add the function
      int idx = int(functions_.size());
      outArg.view(argFunction_).setConstant(idx);
      outDer.view(derFunction_).setConstant(idx);
      functions_.push_back (Function(f, inArg, outArg, inDer, outDer));

      /// Computation order
      std::size_t order = 0;
      computationOrder_.resize(functions_.size());
      Computed_t computed(functions_.size(), false);
      for(std::size_t i = 0; i < functions_.size(); ++i)
        computeOrder(i, order, computed);
      assert(order == functions_.size());

      // Update the free dofs
      outArgs_.addRow(outIdx.first, outIdx.second);
      outArgs_.updateIndexes<true, true, true>();
      inArgs_ = RowBlockIndexes(BlockIndex::difference(BlockIndex::type(0, argSize_), outArgs_.indexes()));
      outDers_.addRow(outDerIdx.first, outDerIdx.second);
      outDers_.updateIndexes<true, true, true>();
      inDers_ = ColBlockIndexes(BlockIndex::difference(BlockIndex::type(0, derSize_), outDers_.indexes()));

      return true;
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

    void ExplicitSolver::computeFunction(const std::size_t& iF, vectorOut_t arg, Computed_t& computed) const
    {
      if (computed[iF]) return;
      const Function& f = functions_[iF];
      // Computed dependencies
      for (std::size_t i = 0; i < f.inArg.indexes().size(); ++i) {
        for (size_type j = 0; j < f.inArg.indexes()[i].second; ++j) {
          if (argFunction_[f.inArg.indexes()[i].first + j] < 0) continue;
          assert((std::size_t)argFunction_[f.inArg.indexes()[i].first + j] < functions_.size());
          computeFunction(argFunction_[f.inArg.indexes()[i].first + j], arg, computed);
        }
      }
      // Compute this function
      // TODO dynamic allocation could be avoided.
      vector_t result(f.f->outputSize());
      f.f->value(result, f.inArg.view(arg).eval());
      f.outArg.view(arg) = result;
      computed[iF] = true;
    }

    void ExplicitSolver::jacobian(matrixOut_t jacobian, vectorIn_t arg) const
    {
      // TODO this could be done only on the complement of inDers_
      jacobian.setZero();
      MatrixOutView_t (jacobian,
          inDers_.nbIndexes(), inDers_.indexes(),
          inDers_.nbIndexes(), inDers_.indexes()).setIdentity();
      // Compute the function jacobians
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        const Function& f = functions_[i];
        f.jacobian.resize(f.f->outputDerivativeSize(), f.f->inputDerivativeSize());
        f.f->jacobian(f.jacobian, f.inArg.view(arg).eval());
      }
      for(std::size_t i = 0; i < functions_.size(); ++i) {
        computeJacobian(computationOrder_[i], jacobian);
      }
    }

    void ExplicitSolver::computeJacobian(const std::size_t& iF, matrixOut_t J) const
    {
      const Function& f = functions_[iF];
      matrix_t Jg (MatrixOutView_t (J,
            f.inDer.nbIndexes(), f.inDer.indexes(),
            inDers_.nbIndexes(), inDers_.indexes()).eval());
      MatrixOutView_t (J,
          f.outDer.nbIndexes(), f.outDer.indexes(),
          inDers_.nbIndexes(), inDers_.indexes()) = f.jacobian * Jg;
    }

    void ExplicitSolver::computeOrder(const std::size_t& iF, std::size_t& iOrder, Computed_t& computed)
    {
      if (computed[iF]) return;
      const Function& f = functions_[iF];
      for (std::size_t i = 0; i < f.inDer.indexes().size(); ++i) {
        for (size_type j = 0; j < f.inDer.indexes()[i].second; ++j) {
          if (derFunction_[f.inDer.indexes()[i].first + j] < 0) continue;
          assert((std::size_t)derFunction_[f.inDer.indexes()[i].first + j] < functions_.size());
          computeOrder(derFunction_[f.inDer.indexes()[i].first + j], iOrder, computed);
        }
      }
      computationOrder_[iOrder] = iF;
      ++iOrder;
      computed[iF] = true;
    }
  } // namespace constraints
} // namespace hpp
