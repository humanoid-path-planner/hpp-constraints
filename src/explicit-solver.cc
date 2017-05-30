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
        const RowBlockIndexes& outArg)
    {
      assert (outArg.indexes().size() == 1 && "Only contiguous function output is supported.");
      const RowBlockIndexes::BlockIndex_t& outIdx = outArg.indexes()[0];

      // Sanity check: is it explicit ?
      for (std::size_t i = 0; i < inArg.indexes().size(); ++i)
        if (BlockIndex::overlap(inArg.indexes()[i], outIdx))
          return false;
      // Check that no other function already computes its outputs.
      if ((outArg.view(dofFunction_).eval().array() >= 0).any())
        return false;
      // Check that it does not insert a loop.
      std::queue<size_type> idxArg;
      append(inArg, idxArg);
      while (!idxArg.empty()) {
        size_type iArg = idxArg.back();
        idxArg.pop();
        if (dofFunction_[iArg] < 0) continue;
        if (iArg >= outIdx.first && iArg < outIdx.first + outIdx.second) return false;
        const Function& func = functions_[dofFunction_[iArg]];
        append(func.inArg, idxArg);
      }

      /// Computation order
      bool isRoot = ( inArg.nbIndexes() == 0
          ||        ( inArg.view(dofFunction_).eval().array() < 0).all()
          );

      int idx = (int)functions_.size();
      outArg.view(dofFunction_).setConstant(idx);
      if (isRoot) functions_.push_front (Function(f, inArg, outArg));
      else        functions_.push_back  (Function(f, inArg, outArg));


      // Update the free dofs
      RowBlockIndexes nFreeDofs;
      const RowBlockIndexes::BlockIndexes_t& idxes = freeDofs_.indexes();
      for (std::size_t i = 0; i < idxes.size(); ++i) {
        RowBlockIndexes::BlockIndexes_t diffs = BlockIndex::difference(idxes[i], outIdx);
        for (std::size_t j = 0; j < diffs.size(); ++j)
          nFreeDofs.addRow(diffs[j].first, diffs[j].second);
      }
      freeDofs_ = nFreeDofs;

      return true;
    }

    void ExplicitSolver::computeFunction(const std::size_t& iF, vectorOut_t arg, Computed_t& computed) const
    {
      if (computed[iF]) return;
      const Function& f = functions_[iF];
      // Computed dependencies
      for (std::size_t i = 0; i < f.inArg.indexes().size(); ++i) {
        for (size_type j = 0; j < f.inArg.indexes()[i].second; ++j) {
          assert(dofFunction_[f.inArg.indexes()[i].first + j] >= 0
              && dofFunction_[f.inArg.indexes()[i].first + j] < functions_.size());
          computeFunction(dofFunction_[f.inArg.indexes()[i].first + j], arg, computed);
        }
      }
      // Compute this function
      // TODO dynamic allocation could be avoided.
      vector_t result(f.f->outputSize());
      f.f->value(result, f.inArg.view(arg).eval());
      f.outArg.view(arg) = result;
      computed[iF] = true;
    }
  } // namespace constraints
} // namespace hpp
