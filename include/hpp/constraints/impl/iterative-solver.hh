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

#ifndef HPP_CONSTRAINTS_IMPL_ITERATIVE_SOLVER_HH
#define HPP_CONSTRAINTS_IMPL_ITERATIVE_SOLVER_HH

namespace hpp {
  namespace constraints {
    namespace lineSearch {
      inline Backtracking::Backtracking () : c (0.001), tau (0.7), smallAlpha (0.2) {}

      template <typename SolverType>
      inline bool Backtracking::operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t u)
      {
        arg_darg.resize(arg.size());

        const value_type slope = computeLocalSlope(solver);
        const value_type t = 2 * c * slope;
        const value_type f_arg_norm2 = solver.residualError();

        if (t > 0) {
          hppDout (error, "The descent direction is not valid: " << t/c);
        } else {
          value_type alpha = 1;

          while (alpha > smallAlpha) {
            darg = alpha * u;
            solver.integration() (arg, darg, arg_darg);
            solver.template computeValue<false> (arg_darg);
            solver.computeError ();
            // Check if we are doing better than the linear approximation with coef
            // multiplied by c < 1
            // t < 0 must hold
            const value_type f_arg_darg_norm2 = solver.residualError();
            if (f_arg_norm2 - f_arg_darg_norm2 >= - alpha * t) {
              arg = arg_darg;
              u = darg;
              return true;
            }
            // Prepare next step
            alpha *= tau;
          }
          hppDout (error, "Could find alpha such that ||f(q)||**2 + "
              << c << " * 2*(f(q)^T * J * dq) is doing worse than "
              "||f(q + alpha * dq)||**2");
        }

        u *= smallAlpha;
        solver.integration() (arg, darg, arg);
        return false;
      }

      template <typename SolverType>
      inline value_type Backtracking::computeLocalSlope(const SolverType& solver) const
      {
        value_type slope = 0;
        for (std::size_t i = 0; i < solver.stacks_.size (); ++i) {
          typename SolverType::Data& d = solver.datas_[i];
          const size_type nrows = d.reducedJ.rows();
          if (df.size() < nrows) df.resize(nrows);
          df.head(nrows).noalias() = d.reducedJ * solver.dqSmall_;
          slope += df.head(nrows).dot(d.error);
        }
        return slope;
      }

      inline FixedSequence::FixedSequence() : alpha (.2), alphaMax (.95), K (.8) {}

      template <typename SolverType>
      inline bool FixedSequence::operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg)
      {
        darg *= alpha;
        alpha = alphaMax - K * (alphaMax - alpha);
        solver.integration() (arg, darg, arg);
        return true;
      }

      inline ErrorNormBased::ErrorNormBased(value_type alphaMin, value_type _a, value_type _b)
          : C (0.5 + alphaMin / 2), K ((1 - alphaMin) / 2), a (_a), b (_b)
      {}

      template <typename SolverType>
      inline bool ErrorNormBased::operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg)
      {
        const value_type r = solver.residualError() / solver.squaredErrorThreshold();
        const value_type alpha = C - K * std::tanh(a * r + b);
        darg *= alpha;
        solver.integration() (arg, darg, arg);
        return true;
      }
    }

    template <typename LineSearchType>
    HierarchicalIterativeSolver::Status HierarchicalIterativeSolver::solve (
        vectorOut_t arg,
        LineSearchType lineSearch) const
    {
      hppDout (info, "before projection: " << arg.transpose ());
      assert (!arg.hasNaN());

      size_type errorDecreased = 3, iter = 0;
      value_type previousSquaredNorm =
	std::numeric_limits<value_type>::infinity();

      // Fill value and Jacobian
      computeValue<true> (arg);
      computeError();

      while (squaredNorm_ > squaredErrorThreshold_ && errorDecreased &&
	     iter < maxIterations_) {

        computeDescentDirection ();
        lineSearch (*this, arg, dq_);

	computeValue<true> (arg);
        computeError ();

	hppDout (info, "squareNorm = " << squaredNorm_);
	--errorDecreased;
	if (squaredNorm_ < previousSquaredNorm) errorDecreased = 3;
	previousSquaredNorm = squaredNorm_;
	++iter;

      }

      hppDout (info, "number of iterations: " << iter);
      if (squaredNorm_ > squaredErrorThreshold_) {
	hppDout (info, "Projection failed.");
        return (!errorDecreased) ? ERROR_INCREASED : MAX_ITERATION_REACHED;
      }
      hppDout (info, "After projection: " << arg.transpose ());
      assert (!arg.hasNaN());
      return SUCCESS;
    }

    template <bool ComputeJac>
    inline void HierarchicalIterativeSolver::computeValue (vectorIn_t arg) const
    {
      for (std::size_t i = 0; i < stacks_.size (); ++i) {
        const DifferentiableFunctionStack& f = stacks_[i];
        Data& d = datas_[i];

        f.value   (d.value, arg);
        if (ComputeJac) f.jacobian(d.jacobian, arg);
        // TODO (*(*it)->comparisonType ()) (v, jacobian);
        d.error = d.value - d.rightHandSide;

        // Copy columns that are not reduced
        if (ComputeJac) reduction_.view (d.jacobian).writeTo(d.reducedJ);
      }
    }

  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_IMPL_HYBRID_SOLVER_HH
