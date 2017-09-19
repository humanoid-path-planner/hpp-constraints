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

#include <pinocchio/multibody/joint/joint.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/liegroup.hh>

namespace se3 {
  using ::hpp::constraints::vectorIn_t;
  using ::hpp::constraints::vectorOut_t;
  using ::hpp::constraints::size_type;

  struct DifferenceStep : public fusion::JointModelVisitor<DifferenceStep >
  {
    typedef boost::fusion::vector<vectorIn_t,
                                  vectorIn_t,
                                  size_type&,
                                  vectorOut_t,
                                  size_type& > ArgsType;

    JOINT_MODEL_VISITOR_INIT(DifferenceStep);

    template<typename JointModel>
    static void algo(const JointModelBase<JointModel>& jmodel,
                    vectorIn_t  q0,
                    vectorIn_t  q1,
                    size_type& rowArg,
                    vectorOut_t result,
                    size_type& rowDer)
    {
      ::hpp::pinocchio::LieGroupTpl::template operation<JointModel>::type ::difference (
          q0.segment(rowArg, jmodel.nq()),
          q1.segment(rowArg, jmodel.nq()),
          result.segment(rowDer, jmodel.nv()));
      rowArg += jmodel.nq();
      rowDer += jmodel.nv();
    }
  };

  template<>
  void DifferenceStep::algo (const JointModelBase<JointModelComposite>& jmodel,
                    vectorIn_t  q0,
                    vectorIn_t  q1,
                    size_type& rowArg,
                    vectorOut_t result,
                    size_type& rowDer)
  {
    ::se3::details::Dispatch< DifferenceStep >::run(jmodel.derived(),
        ArgsType(q0, q1, rowArg, result, rowDer));
  }
}


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

    void difference (const DevicePtr_t& robot,
        const Eigen::BlockIndex<size_type>::vector_t& indexes,
        vectorIn_t arg0,
        vectorIn_t arg1,
        vectorOut_t result)
    {
      typedef typename se3::DifferenceStep DiffStep;
      const se3::Model& model = robot->model();
      std::size_t iJoint = 1;

      size_type rowArg = 0, rowDer = 0;
      for (std::size_t i = 0; i < indexes.size(); ++i) {
        const Eigen::BlockIndex<size_type>::type& interval = indexes[i];
        size_type j = 0;
        while (j < interval.second) {
          size_type iArg = interval.first + j;

          if (iArg >= model.nq) { // Extra dofs, assume vector space
            // TODO this could be optimized for cases where there are many
            // extra dofs.
            result[rowDer] = arg1[rowArg] - arg0[rowArg];
            ++rowArg;
            ++rowDer;
            continue;
          }
          while (model.joints[iJoint].idx_q() != iArg) {
            ++iJoint;
            if (iJoint >= model.joints.size()) throw std::runtime_error("Joint index out of bounds");
          }

          const se3::JointModel& jmodel = model.joints[iJoint];

          assert (jmodel.idx_q() >= interval.first);
          assert (jmodel.nq()    <= interval.second);

          // Compute difference
          DiffStep::ArgsType args (arg0, arg1, rowArg, result, rowDer);
          DiffStep::run(model.joints[iJoint], args);

          ++iJoint;
          j += jmodel.nq();
        }
      }
    }

    Eigen::ColBlockIndexes ExplicitSolver::activeParameters () const
    {
      BlockIndex::vector_t biv;
      for (std::size_t i = 0; i < functions_.size (); ++i)
        biv.insert(biv.end(), functions_[i].inArg.indexes().begin(),
                              functions_[i].inArg.indexes().end());
      ColBlockIndexes cbi (biv);
      cbi.updateIndexes<true, true, true>();
      return cbi;
    }

    Eigen::ColBlockIndexes ExplicitSolver::activeDerivativeParameters () const
    {
      BlockIndex::vector_t biv;
      for (std::size_t i = 0; i < functions_.size (); ++i)
        biv.insert(biv.end(), functions_[i].inDer.indexes().begin(),
                              functions_[i].inDer.indexes().end());
      ColBlockIndexes cbi (biv);
      cbi.updateIndexes<true, true, true>();
      return cbi;
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
      assert(error.size() == outDers_.nbIndexes());
      arg_ = arg;
      solve (arg_);
      difference_ (arg, arg_, diff_);
      outDers_.rview(diff_).writeTo(error);
      hppDout (info, "Squared error norm is " << error.squaredNorm());
      // TODO: this threshold a very bad way of solving a numerical issue.
      // in hpp-core, a numerical constraint may have both a explicit and an
      // implicit formulation.
      // 1. The implicit formulation never gives an exact solution so we must
      //    allow a threshold.
      // 2. The explicit function may give slightly different results so we must
      //    allow for a greater threshold. Note the factor 10 is very likely
      //    over-estimated.
      return error.squaredNorm() < squaredErrorThreshold_;
    }

    bool ExplicitSolver::isSatisfied (vectorIn_t arg) const
    {
      diffSmall_.resize(outDers_.nbIndexes());
      return isSatisfied (arg, diffSmall_);
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
      if ((outArg.rview(argFunction_).eval().array() >= 0).any())
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
      outArg.lview(argFunction_).setConstant(idx);
      outDer.lview(derFunction_).setConstant(idx);
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

    void ExplicitSolver::computeFunction(const std::size_t& iF, vectorOut_t arg) const
    {
      const Function& f = functions_[iF];
      // Compute this function
      f.f->value(f.value, f.inArg.rview(arg).eval());
      f.outArg.lview(arg) = f.value.vector ();
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
        f.f->jacobian(f.jacobian, f.inArg.rview(arg).eval());
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
