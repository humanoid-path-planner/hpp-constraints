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

#include <hpp/constraints/matrix-view.hh>
#include <hpp/constraints/differentiable-function-stack.hh>

namespace hpp {
  namespace constraints {
    namespace lineSearch {
      /// Implements the backtracking line search algorithm
      /// See https://en.wikipedia.org/wiki/Backtracking_line_search
      struct Backtracking {
        Backtracking ();

        template <typename SolverType>
        inline bool operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg);

        template <typename SolverType>
        inline value_type computeLocalSlope(const SolverType& solver) const;

        const value_type c, tau, smallAlpha; // 0.8 ^ 7 = 0.209, 0.8 ^ 8 = 0.1677
        mutable vector_t arg_darg, df, darg;
      };

      /// The step size is computed using the recursion:
      /// \f[ alpha <- alpha - K * (alphaMax - alpha) \f]
      struct FixedSequence {
        FixedSequence();

        template <typename SolverType>
        inline bool operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg);

        value_type alpha;
        const value_type alphaMax, K;
      };

      /// The step size is computed using the formula
      /// \f[ const value_type alpha = C - K * std::tanh(a * r + b) \f]
      struct ErrorNormBased {
        enum { N = 4, M = 8 };

        ErrorNormBased(value_type alphaMin = 0.2,
            value_type _a = 4 / (std::pow(10, (int)M) - std::pow(10, (int)N)),
            value_type _b = 2 - 4 / (1 - std::pow (10, (int)(N - M))));

        template <typename SolverType>
        inline bool operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg);

        const value_type C, K, a, b;
      };
    }

    class HPP_CONSTRAINTS_DLLAPI HierarchicalIterativeSolver
    {
      public:
        typedef Eigen::ColBlockIndexes Reduction_t;
        typedef lineSearch::FixedSequence DefaultLineSearch;

        enum Status {
          ERROR_INCREASED,
          MAX_ITERATION_REACHED,
          SUCCESS
        };
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

        /// Set the velocity variable that must be changed.
        /// The other variables will be left unchanged by the iterative
        /// algorithm.
        void reduction (const intervals_t intervals)
        {
          reduction_ = Reduction_t();
          for (std::size_t i = 0; i < intervals.size(); ++i)
            reduction_.addCol(intervals[i].first, intervals[i].second);
          reduction_.updateIndexes<true, true, true>();
          update ();
        }

        /// Set the velocity variable that must be changed.
        /// The other variables will be left unchanged by the iterative
        /// algorithm.
        void reduction (const Reduction_t& reduction)
        {
          reduction_ = reduction;
          update ();
        }

        template <typename LineSearchType>
        Status solve (vectorOut_t arg, LineSearchType ls = LineSearchType()) const;

        inline Status solve (vectorOut_t arg) const
        {
          return solve (arg, DefaultLineSearch());
        }

        /// Set the integration function
        void integration (const Integration_t& integrate)
        {
          integrate_ = integrate;
        }

        /// Get the integration function
        const Integration_t& integration () const
        {
          return integrate_;
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
        /// Get error threshold
        value_type errorThreshold () const
        {
          return sqrt (squaredErrorThreshold_);
        }
        /// Get error threshold
        value_type squaredErrorThreshold () const
        {
          return squaredErrorThreshold_;
        }

        value_type residualError() const
        {
          return squaredNorm_;
        }

        void lastIsOptional (bool optional)
        {
          lastIsOptional_ = optional;
        }

        bool lastIsOptional () const
        {
          return lastIsOptional_;
        }

        /// Compute a right hand side using the input arg.
        /// This does not set the right hand side.
        /// To set the right hand side using this function, one must call
        /// rightHandSide (rightHandSideFromArgument (arg))
        vector_t rightHandSideFromInput (vectorIn_t arg) const;

        /// Set the level set parameter.
        /// \param rhs the level set parameter.
        void rightHandSide (const vector_t& rhs);

        /// Get the level set parameter.
        /// \return the parameter.
        vector_t rightHandSide () const;

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

        template <bool ComputeJac> void computeValue (vectorIn_t arg) const;
        void computeError () const;
        void computeDescentDirection () const;
        void expandDqSmall () const;

        value_type squaredErrorThreshold_;
        size_type maxIterations_;

        std::vector<DifferentiableFunctionStack> stacks_;
        size_type dimension_;
        bool lastIsOptional_;
        Reduction_t reduction_;
        Integration_t integrate_;

        mutable vector_t dq_, dqSmall_;
        mutable matrix_t projector_;
        mutable value_type squaredNorm_;
        mutable std::vector<Data> datas_;

        mutable ::hpp::statistics::SuccessStatistics statistics_;

        friend struct lineSearch::Backtracking;
    }; // class IterativeSolver
  } // namespace constraints
} // namespace hpp

#include <hpp/constraints/impl/iterative-solver.hh>

#endif // HPP_CONSTRAINTS_ITERATIVE_SOLVER_HH
