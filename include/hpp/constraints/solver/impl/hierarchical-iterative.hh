// Copyright (c) 2017, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#ifndef HPP_CONSTRAINTS_SOLVER_IMPL_HIERARCHICAL_ITERATIVE_HH
#define HPP_CONSTRAINTS_SOLVER_IMPL_HIERARCHICAL_ITERATIVE_HH

#include <hpp/util/debug.hh>

#include <hpp/constraints/config.hh>
#include <hpp/constraints/svd.hh>

namespace hpp {
  namespace constraints {
    namespace solver {
    namespace lineSearch {
      template <typename SolverType>
      inline bool Constant::operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg)
      {
        solver.integrate (arg, darg, arg);
        return true;
      }

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
          /* TODO add a threshold to avoid too large steps.
          const value_type u2 = u.squaredNorm();
          if (u2 > 1.) alpha = 1. / std::sqrt(u2);
          */

          while (alpha > smallAlpha) {
            darg = alpha * u;
            solver.integrate (arg, darg, arg_darg);
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
        solver.integrate (arg, u, arg);
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
          slope += df.head(nrows).dot(d.activeRowsOfJ.keepRows().rview(d.error).eval());
        }
        return slope;
      }

      template <typename SolverType>
      inline bool FixedSequence::operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg)
      {
        darg *= alpha;
        alpha = alphaMax - K * (alphaMax - alpha);
        solver.integrate (arg, darg, arg);
        return true;
      }

      template <typename SolverType>
      inline bool ErrorNormBased::operator() (const SolverType& solver, vectorOut_t arg, vectorOut_t darg)
      {
        const value_type r = solver.residualError() / solver.squaredErrorThreshold();
        const value_type alpha = C - K * std::tanh(a * r + b);
        darg *= alpha;
        solver.integrate (arg, darg, arg);
        return true;
      }
    }

    template <typename LineSearchType>
    inline solver::HierarchicalIterative::Status solver::HierarchicalIterative::solve (
        vectorOut_t arg,
        LineSearchType lineSearch) const
    {
      hppDout (info, "before projection: " << arg.transpose ());
      assert (!arg.hasNaN());

      size_type errorDecreased = 3, iter = 0;
      value_type previousSquaredNorm =
	std::numeric_limits<value_type>::infinity();
      static const value_type dqMinSquaredNorm = Eigen::NumTraits<value_type>::dummy_precision();

      // Fill value and Jacobian
      computeValue<true> (arg);
      computeError();

      if (squaredNorm_ > squaredErrorThreshold_
          && reducedDimension_ == 0) return INFEASIBLE;

      Status status;
      while (squaredNorm_ > squaredErrorThreshold_ && errorDecreased &&
	     iter < maxIterations_) {

        computeSaturation(arg);
        computeDescentDirection ();
        if (dq_.squaredNorm () < dqMinSquaredNorm) {
          // TODO INFEASIBLE means that we have reached a local minima.
          // The problem may still be feasible from a different starting point.
          status = INFEASIBLE;
          break;
        }
        lineSearch (*this, arg, dq_);

	computeValue<true> (arg);
        computeError ();

	hppDout (info, "squareNorm = " << squaredNorm_);
	--errorDecreased;
	if (squaredNorm_ < previousSquaredNorm)
          errorDecreased = 3;
        else
          status = ERROR_INCREASED;
	previousSquaredNorm = squaredNorm_;
	++iter;

      }

      hppDout (info, "number of iterations: " << iter);
      if (squaredNorm_ > squaredErrorThreshold_) {
	hppDout (info, "Projection failed.");
        return (iter >= maxIterations_) ? MAX_ITERATION_REACHED : status;
      }
      hppDout (info, "After projection: " << arg.transpose ());
      assert (!arg.hasNaN());
      return SUCCESS;
    }
    } // namespace solver
  } // namespace constraints
} // namespace hpp

#endif //HPP_CONSTRAINTS_SOLVER_IMPL_HIERARCHICAL_ITERATIVE_HH
