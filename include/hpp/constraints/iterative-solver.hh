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
      /// No line search. Use \f$\alpha_i = 1\f$
      struct Constant {
        template <typename SolverType>
        bool operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg);
      };

      /// Implements the backtracking line search algorithm.
      /// See https://en.wikipedia.org/wiki/Backtracking_line_search.
      struct Backtracking {
        Backtracking ();

        template <typename SolverType>
        bool operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg);

        template <typename SolverType>
        inline value_type computeLocalSlope(const SolverType& solver) const;

        value_type c, tau, smallAlpha; // 0.8 ^ 7 = 0.209, 0.8 ^ 8 = 0.1677
        mutable vector_t arg_darg, df, darg;
      };

      /// The step size is computed using the recursion
      /// \f$ \alpha_{i+1} = \alpha - K \times (\alpha_{max} - \alpha_i) \f$
      /// where \f$K\f$ and \f$\alpha_{max}\f$ are some constant values.
      struct FixedSequence {
        FixedSequence();

        template <typename SolverType>
        bool operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg);

        value_type alpha;
        value_type alphaMax, K;
      };

      /// The step size is computed using the formula
      /// \f$ \alpha_{i} = C - K \times \text{tanh}(a \frac{\|f(\mathbf{q}_i)\|}{\epsilon^2} + b) \f$, where
      /// \li \f$\epsilon\f$ is the error threshold:
      /// if \f$\|f(\mathbf{q}_i)\|<\epsilon\f$, \f$\mathbf{q}_i\f$ is
      /// considered to satisfy the constraint.
      struct ErrorNormBased {
        ErrorNormBased(value_type alphaMin, value_type _a, value_type _b);
        ErrorNormBased(value_type alphaMin = 0.2);

        template <typename SolverType>
        bool operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg);

        value_type C, K, a, b;
      };
    }

    /// Solve a system of non-linear equations on a robot configuration
    ///
    /// The non-linear system of equations is built by adding equations with
    /// method HierarchicalIterativeSolver::add.
    ///
    /// Note that a hierarchy between the equations can be
    /// provided. In this case, the solver will try to solve the
    /// highest priority equations first, and then to solve the lower priority
    /// equations. Note that priorities are in decreasing order: 0 has higher
    /// priority than 1.
    ///
    /// The algorithm used is a Newton-Raphson like algorithm that works as
    /// follows: let \f$f (\mathbf{q}) = 0\f$ be the system of equations where
    /// \f$f\f$ is a \f$C^1\f$ mapping from the robot configuration space to
    /// a Lie group space \f$\mathcal{L}\f$.
    ///
    /// Starting from initial guess \f$\mathbf{q}_0\f$, the method
    /// HierarchicalIterativeSolver::solve builds a sequence of configurations
    /// \f$\mathbf{q}_i\f$ as follows:
    /// \f{eqnarray*}
    /// \mathbf{q}_{i+1} = \mathbf{q}_i -
    ///    \alpha_i \frac{\partial f}{\partial \mathbf{q}}(\mathbf{q}_i)^{+}
    ///    f (\mathbf{q}_i)
    /// \f}
    /// where
    /// \li \f$\frac{\partial f}{\partial \mathbf{q}}(\mathbf{q}_i)^{+}\f$ is
    ///     the Moore-Penrose pseudo-inverse of the system Jacobian,
    /// \li \f$\alpha_i\f$ is a sequence of real numbers depending on the
    ///     line search strategy. Possible line-search strategies are
    ///     lineSearch::Constant, lineSearch::Backtracking,
    ///     lineSearch::FixedSequence, lineSearch::ErrorNormBased.
    /// until
    /// \li the residual \f$\|f(\mathbf{q})\|\f$ is below an error threshold, or
    /// \li the maximal number of iterations has been reached.
    ///
    /// The error threshold can be accessed by methods
    /// HierarchicalIterativeSolver::errorThreshold. The maximal number of
    /// iterations can be accessed by methods
    /// HierarchicalIterativeSolver::maxIterations.
    ///
    /// \note Lie group
    ///
    /// The unknowns \f$\mathbf{q}\f$ may take values in a more general set
    /// than the configuration space of a robot. This set should be a Cartesian
    /// product of Lie groups. In this case, the user can provide a method that
    /// computes the exponential map of a tangent vector.
    /// \sa HierarchicalIterativeSolver::Integration_t and
    /// HierarchicalIterativeSolver::integration.
    ///
    /// \note Saturation
    ///
    /// To prevent configuration variables to get out of joint limits during
    /// Newton Raphson iterations, the user may provide a method of type
    /// HierarchicalIterativeSolver::Saturation_t using setter and getter
    /// HierarchicalIterativeSolver::saturation.
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
        /// This function checks which degrees of freedom are saturated.
        ///
        /// \param result a configuration
        ///
        /// For each degree of freedom, saturation is set to
        /// \li -1 if the lower bound is reached,
        /// \li  1 if the upper bound is reached,
        /// \li  0 otherwise.
        typedef boost::function<bool (vectorIn_t result, Eigen::VectorXi& saturation)> Saturation_t;

        HierarchicalIterativeSolver (const std::size_t& argSize, const std::size_t derSize);

        virtual ~HierarchicalIterativeSolver () {}

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

        /// Set the saturation function
        void saturation (const Saturation_t& saturate)
        {
          saturate_ = saturate;
        }

        /// Get the saturation function
        const Saturation_t& saturation () const
        {
          return saturate_;
        }

        /// \}

        /// \name Problem resolution
        /// \{

        /// Solve the system of non linear equations
        ///
        /// \param arg initial guess,
        /// \param ls line search method used.
        ///
        /// Use Newton Rhapson like iterative method until the error is below
        /// the threshold, or until the maximal number of iterations has been
        /// reached.
        ///
        /// \note Explicit constraints are expressed in their implicit
        ///       form: \f$\mathbf{q}_2 = f (\mathbf{q}_1)\f$ is replaced by
        ///       \f$\mathbf{q}_2 - f (\mathbf{q}_1) = 0\f$.
        template <typename LineSearchType>
        Status solve (vectorOut_t arg, LineSearchType ls = LineSearchType()) const;

        /// Solve the system of non linear equations
        ///
        /// \param arg initial guess,
        ///
        /// Use Newton Rhapson like iterative method until the error is below
        /// the threshold, or until the maximal number of iterations has been
        /// reached. Use the default line search method (fixed sequence of
        /// \f$\alpha_i\f$).
        ///
        /// \note Explicit constraints are expressed in their implicit
        ///       form.
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

        /// Dimension of the problem after removing the rows of the jacobian
        /// which do not influence the error (only zeros along those lines).
        const size_type& reducedDimension () const
        {
          return reducedDimension_;
        }

        /// Configuration parameters involved in the constraint resolution.
        ArrayXb activeParameters () const;

        /// Velocity parameters involved in the constraint resolution.
        ArrayXb activeDerivativeParameters () const;

        /// \}

        /// Returns the squared norm of the error vector
        value_type residualError() const
        {
          return squaredNorm_;
        }

        /// Returns the error vector
        void residualError(vectorOut_t error) const;

        /// \name Right hand side accessors
        /// \{

        /// Compute a right hand side using the input arg.
        vector_t rightHandSideFromInput (vectorIn_t arg);

        bool rightHandSideFromInput (const DifferentiableFunctionPtr_t& f, vectorIn_t arg);

        bool rightHandSide (const DifferentiableFunctionPtr_t& f, vectorIn_t rhs);

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
        void computeSaturation (vectorIn_t arg) const;
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

        virtual void integrate(vectorIn_t from, vectorIn_t velocity, vectorOut_t result) const
        {
          integrate_ (from, velocity, result);
        }

        /// \}

        virtual std::ostream& print (std::ostream& os) const;

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
        void saturate (vectorOut_t arg) const;


        value_type squaredErrorThreshold_, inequalityThreshold_;
        size_type maxIterations_;

        std::vector<DifferentiableFunctionStack> stacks_;
        size_type argSize_, derSize_;
        size_type dimension_, reducedDimension_;
        bool lastIsOptional_;
        Reduction_t reduction_;
        Integration_t integrate_;
        Saturation_t saturate_;
        /// The smallest non-zero singular value
        mutable value_type sigma_;

        mutable vector_t dq_, dqSmall_;
        mutable matrix_t projector_, reducedJ_;
        mutable Eigen::VectorXi saturation_, reducedSaturation_;
        mutable ArrayXb tmpSat_;
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
