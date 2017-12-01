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

    /// \addtogroup solvers
    /// \{
    namespace lineSearch {
      /// No line search. Use alpha = 1
      struct Constant {
        template <typename SolverType>
        bool operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg);
      };

      /// Implements the backtracking line search algorithm
      /// See https://en.wikipedia.org/wiki/Backtracking_line_search
      struct Backtracking {
        Backtracking ();

        template <typename SolverType>
        bool operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg);

        template <typename SolverType>
        inline value_type computeLocalSlope(const SolverType& solver) const;

        value_type c, tau, smallAlpha; // 0.8 ^ 7 = 0.209, 0.8 ^ 8 = 0.1677
        mutable vector_t arg_darg, df, darg;
      };

      /// The step size is computed using the recursion:
      /// \f[ alpha <- alpha - K * (alphaMax - alpha) \f]
      struct FixedSequence {
        FixedSequence();

        template <typename SolverType>
        bool operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg);

        value_type alpha;
        value_type alphaMax, K;
      };

      /// The step size is computed using the formula
      /// \f[ const value_type alpha = C - K * std::tanh(a * r + b) \f]
      struct ErrorNormBased {
        ErrorNormBased(value_type alphaMin, value_type _a, value_type _b);
        ErrorNormBased(value_type alphaMin = 0.2);

        template <typename SolverType>
        bool operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg);

        value_type C, K, a, b;
      };
    }

    class HPP_CONSTRAINTS_DLLAPI HierarchicalIterativeSolver
    {
      public:
        typedef Eigen::ColBlockIndices Reduction_t;
        typedef lineSearch::FixedSequence DefaultLineSearch;

        enum Status {
          ERROR_INCREASED,
          MAX_ITERATION_REACHED,
          INFEASIBLE,
          SUCCESS
        };
        /// This function integrates velocity during unit time, from argument.
        /// It should be robust to cases where from and result points to the
        /// same vector in memory (aliasing)
        typedef boost::function<void (vectorIn_t from, vectorIn_t velocity, vectorOut_t result)> Integration_t;

        HierarchicalIterativeSolver (const std::size_t& argSize, const std::size_t derSize);

        /// \name Problem definition
        /// \{

        /// Add an equality constraint to a priority
        void add (const DifferentiableFunctionPtr_t& f, const std::size_t& priority)
        {
          add (f, priority, ComparisonTypes_t(f->outputSize(), EqualToZero));
        }

        /// Add a constraint \f{ f(q) comp 0 \f} to a priority
        void add (const DifferentiableFunctionPtr_t& f, const std::size_t& priority, const ComparisonTypes_t& comp);

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

        /// \}

        /// \name Problem resolution
        /// \{

        template <typename LineSearchType>
        Status solve (vectorOut_t arg, LineSearchType ls = LineSearchType()) const;

        inline Status solve (vectorOut_t arg) const
        {
          return solve (arg, DefaultLineSearch());
        }

        bool isSatisfied (vectorIn_t arg) const
        {
          computeValue<false>(arg);
          computeError();
          return squaredNorm_ < squaredErrorThreshold_;
        }

        /// Returns the lowest singular value.
        /// If the jacobian has maximum rank r, then it corresponds to r-th
        /// greatest singular value. This value is zero when the jacobian is
        /// singular.
        const value_type& sigma () const
        {
          return sigma_;
        }

        /// \}

        /// \name Parameters
        /// \{

        /// Set the velocity variable that must be changed.
        /// The other variables will be left unchanged by the iterative
        /// algorithm.
        void reduction (const segments_t intervals)
        {
          reduction_ = Reduction_t();
          for (std::size_t i = 0; i < intervals.size(); ++i)
            reduction_.addCol(intervals[i].first, intervals[i].second);
          reduction_.updateIndices<true, true, true>();
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

        /// Get the inequality threshold
        value_type inequalityThreshold () const
        {
          return inequalityThreshold_;
        }
        /// set the inequality threshold
        void inequalityThreshold (const value_type& it)
        {
          inequalityThreshold_ = it;
        }

        void lastIsOptional (bool optional)
        {
          lastIsOptional_ = optional;
        }

        bool lastIsOptional () const
        {
          return lastIsOptional_;
        }

        /// \}

        /// \name Stack
        /// \{

        const DifferentiableFunctionStack& stack(const std::size_t priority)
        {
          assert(priority < stacks_.size());
          return stacks_[priority];
        }

        std::size_t numberStacks() const
        {
          return stacks_.size();
        }

        const size_type& dimension () const
        {
          return dimension_;
        }

        /// Configuration parameters involved in the constraint resolution.
        bool_array_t activeParameters () const;

        /// Velocity parameters involved in the constraint resolution.
        bool_array_t activeDerivativeParameters () const;

        /// \}

        /// Returns the squared norm of the error vector
        value_type residualError() const
        {
          return squaredNorm_;
        }

        /// Returns the error vector
        void residualError(vectorOut_t error) const
        {
          size_type row = 0;
          for (std::size_t i = 0; i < datas_.size(); ++i) {
            const Data& d = datas_[i]; 
            error.segment(row, d.error.size()) = d.error;
            row += d.error.size();
          }
        }

        /// \name Right hand side accessors
        /// \{

        /// Compute a right hand side using the input arg.
        vector_t rightHandSideFromInput (vectorIn_t arg);

        bool rightHandSideFromInput (const DifferentiableFunctionPtr_t& f, vectorIn_t arg) const;

        bool rightHandSide (const DifferentiableFunctionPtr_t& f, vectorIn_t rhs) const;

        /// Set the level set parameter.
        /// \param rhs the level set parameter.
        void rightHandSide (vectorIn_t rhs);

        /// Get the level set parameter.
        /// \return the parameter.
        vector_t rightHandSide () const;

        /// Get size of the level set parameter.
        size_type rightHandSideSize () const;

        /// \}

        /// \name Access to internal datas
        /// You should know what you do when you call these functions
        /// \{

        /// Compute the value of each level, and the jacobian if ComputeJac is true.
        template <bool ComputeJac> void computeValue (vectorIn_t arg) const;
        void getValue (vectorOut_t v) const;
        void getReducedJacobian (matrixOut_t J) const;
        /// If lastIsOptional() is true, then the last level is ignored.
        /// \warning computeValue must have been called first.
        void computeError () const;

        /// Accessor to the last step done
        const vector_t& lastStep () const
        {
          return dq_;
        }
        /// \}

      protected:
        typedef Eigen::JacobiSVD <matrix_t> SVD_t;

        struct Data {
          /// \cond
          EIGEN_MAKE_ALIGNED_OPERATOR_NEW
          /// \endcond
          LiegroupElement output, rightHandSide;
          vector_t error;
          matrix_t jacobian, reducedJ;

          SVD_t svd;
          matrix_t PK;

          mutable size_type maxRank;

          ComparisonTypes_t comparison;
          std::vector<std::size_t> inequalityIndices;
          Eigen::RowBlockIndices equalityIndices;
          Eigen::MatrixBlocks<false,false> activeRowsOfJ;
        };

        /// Allocate datas and update sizes of the problem
        /// Should be called whenever the stack is modified.
        void update ();

        /// Compute which rows of the jacobian of stack_[iStack]
        /// are not zero, using the activeDerivativeParameters of the functions.
        /// The result is stored in datas_[i].activeRowsOfJ
        virtual void computeActiveRowsOfJ (std::size_t iStack);

        /// Compute a SVD decomposition of each level and find the best descent
        /// direction at the first order.
        /// Linearization of the system of equations
        /// rhs - v_{i} = J (q_i) (dq_{i+1} - q_{i})
        /// q_{i+1} - q_{i} = J(q_i)^{+} ( rhs - v_{i} )
        /// dq = J(q_i)^{+} ( rhs - v_{i} )
        /// \warning computeValue<true> must have been called first.
        void computeDescentDirection () const;
        void expandDqSmall () const;
        void integrate(vectorIn_t from, vectorIn_t velocity, vectorOut_t result) const
        {
          integrate_ (from, velocity, result);
        }

        value_type squaredErrorThreshold_, inequalityThreshold_;
        size_type maxIterations_;

        std::vector<DifferentiableFunctionStack> stacks_;
        size_type argSize_, derSize_;
        size_type dimension_;
        bool lastIsOptional_;
        Reduction_t reduction_;
        Integration_t integrate_;
        /// The smallest non-zero singular value
        mutable value_type sigma_;

        mutable vector_t dq_, dqSmall_;
        mutable matrix_t projector_, reducedJ_;
        mutable value_type squaredNorm_;
        mutable std::vector<Data> datas_;
        mutable SVD_t svd_;

        mutable ::hpp::statistics::SuccessStatistics statistics_;

        friend struct lineSearch::Backtracking;
    }; // class IterativeSolver
    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_ITERATIVE_SOLVER_HH
