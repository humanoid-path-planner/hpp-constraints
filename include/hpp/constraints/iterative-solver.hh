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

#ifndef HPP_CONSTRAINTS_ITERATIVE_SOLVER_HH
#define HPP_CONSTRAINTS_ITERATIVE_SOLVER_HH

#include <boost/function.hpp>

#include <hpp/statistics/success-bin.hh>

#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/config.hh>

#include <hpp/constraints/solver.hh>
#include <hpp/constraints/differentiable-function-stack.hh>

namespace hpp {
  namespace constraints {
    class HPP_CONSTRAINTS_DLLAPI HierarchicalIterativeSolver : public Solver
    {
      public:
        /// This function integrates velocity during unit time, from argument.
        /// It should be robust to cases where from and result points to the
        /// same vector in memory (aliasing)
        typedef boost::function<void (vectorIn_t from, vectorIn_t velocity, vectorOut_t result)> Integration_t;

        HierarchicalIterativeSolver ();

        DifferentiableFunctionStack& stack(const std::size_t priority)
        {
          assert(priority < stacks_.size());
          return stacks_[priority];
        }

        std::size_t numberStacks() const
        {
          return stacks_.size();
        }

        void addStack ()
        {
          stacks_.push_back(DifferentiableFunctionStack());
          datas_.push_back(Data());
        }

        void update ();

        void reduction (const intervals_t intervals)
        {
          reduction_ = intervals;
          update ();
        }

        bool solve (vectorOut_t arg) const;

        /// Set the integration function
        void integration (const Integration_t& integrate)
        {
          integrate_ = integrate;
        }

        /// Set maximal number of iterations
        void maxIterations (size_type iterations)
        {
          maxIterations_ = iterations;
        }
        /// Get maximal number of iterations in config projector
        size_type maxIterations () const
        {
          return maxIterations_;
        }

        /// Set error threshold
        void errorThreshold (const value_type& threshold)
        {
          squaredErrorThreshold_ = threshold * threshold;
        }
        /// Get errorimal number of threshold in config projector
        value_type errorThreshold () const
        {
          return sqrt (squaredErrorThreshold_);
        }

        value_type residualError() const
        {
          return squaredNorm_;
        }

      protected:
        typedef Eigen::JacobiSVD <matrix_t> SVD_t;

        struct Data {
          /// \cond
          EIGEN_MAKE_ALIGNED_OPERATOR_NEW
          /// \endcond
          vector_t value, rightHandSide, error;
          matrix_t jacobian, reducedJ;

          SVD_t svd;
          matrix_t PK;
        };

        void computeValueAndJacobian (vectorIn_t arg) const;
        void computeValueAndJacobian (vectorIn_t arg, const std::size_t priority) const;
        void computeError () const;
        void computeIncrement (const value_type& alpha) const;
        void expandDqSmall () const;

        value_type squaredErrorThreshold_;
        size_type maxIterations_;

        std::vector<DifferentiableFunctionStack> stacks_;
        bool lastIsOptional_;
        intervals_t reduction_;
        Integration_t integrate_;

        mutable vector_t dq_, dqSmall_;
        mutable matrix_t projector_;
        mutable value_type squaredNorm_;
        mutable std::vector<Data> datas_;

        mutable ::hpp::statistics::SuccessStatistics statistics_;
    }; // class IterativeSolver
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_ITERATIVE_SOLVER_HH
