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

#ifndef HPP_CONSTRAINTS_SOLVER_IMPL_BY_SUBSTITUTION_HH
#define HPP_CONSTRAINTS_SOLVER_IMPL_BY_SUBSTITUTION_HH

namespace hpp {
  namespace constraints {
    namespace solver {
    typedef std::numeric_limits<value_type> numeric_limits;
    typedef Eigen::NumTraits<value_type> NumTraits;

    template <typename LineSearchType>
    inline HierarchicalIterative::Status BySubstitution::impl_solve (
        vectorOut_t arg,
        bool _optimize,
        LineSearchType lineSearch) const
    {
      bool optimize = _optimize && lastIsOptional_;
      assert (!arg.hasNaN());

      explicit_.solve(arg);

      size_type errorDecreased = 3, iter = 0;
      value_type previousSquaredNorm;
      static const value_type dqMinSquaredNorm = NumTraits::dummy_precision();
      value_type initSquaredNorm = 0;

      // Variables for optimization only
      value_type previousCost = numeric_limits::infinity();
      value_type scaling = 1.;
      bool onlyLineSearch = false;
      vector_t qopt;

      // Fill value and Jacobian
      computeValue<true> (arg);
      computeError();

      bool errorWasBelowThr = (squaredNorm_ < squaredErrorThreshold_);
      vector_t initArg;
      if (errorWasBelowThr) {
        initArg = arg;
        if (!optimize) iter = std::max (maxIterations_,size_type(2)) - 2;
        initSquaredNorm = squaredNorm_;
      }

      bool errorIsAboveThr = (squaredNorm_ > .25 * squaredErrorThreshold_);
      if (errorIsAboveThr && reducedDimension_ == 0) return INFEASIBLE;
      if (optimize && !errorIsAboveThr) qopt = arg;

      Status status = SUCCESS;
      while ( (optimize || (errorIsAboveThr && errorDecreased))) {
        // 1. Maximum iterations
        if (iter >= maxIterations_) {
          status = MAX_ITERATION_REACHED;
          break;
        }
        status = SUCCESS;

        // 2. Compute step
        // onlyLineSearch is true when we only reduced the scaling.
        if (!onlyLineSearch) {
          previousSquaredNorm = squaredNorm_;
          // Update the jacobian using the jacobian of the explicit system.
          updateJacobian(arg);
          computeSaturation(arg);
          computeDescentDirection ();
        }
        // Apply scaling to avoid too large steps.
        if (optimize) dq_ *= scaling;
        if (dq_.squaredNorm () < dqMinSquaredNorm) {
          // We assume that the algorithm reached a local minima.
          status = INFEASIBLE;
          break;
        }
        // 3. Apply line search algorithm for the computed step
        lineSearch (*this, arg, dq_);
        explicit_.solve(arg);

        // 4. Evaluate the error at the new point.
	computeValue<true> (arg);
        computeError ();

	--errorDecreased;
	if (squaredNorm_ < previousSquaredNorm)
          errorDecreased = 3;
        else
          status = ERROR_INCREASED;

        errorIsAboveThr = (squaredNorm_ > .25 * squaredErrorThreshold_);
        // 5. In case of optimization,
        // - if the constraints is satisfied and the cost decreased, increase
        //   the scaling (amount of confidence in the linear approximation)
        // - if the constraints is not satisfied, decrease the scaling and
        //   and cancel this step.
        if (optimize) {
          if (!errorIsAboveThr) {
            value_type cost = datas_.back().error.squaredNorm();
            if (cost < previousCost) {
              qopt = arg;
              previousCost = cost;
              if (scaling < 0.5) scaling *= 2;
            }
            onlyLineSearch = false;
          } else {
            dq_ /= scaling;
            scaling *= 0.5;
            arg = qopt;
            onlyLineSearch = true;
          }
        }

	++iter;
      }

      if (!optimize && errorWasBelowThr) {
        if (squaredNorm_ > initSquaredNorm) {
          arg = initArg;
        }
        return SUCCESS;
      }
      // If optimizing, qopt is the visited configuration that satisfies the
      // constraints and has lowest cost.
      if (optimize && qopt.size() > 0) arg = qopt;

      assert (!arg.hasNaN());
      return status;
    }
    } // namespace solver
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_SOLVER_IMPL_BY_SUBSTITUTION_HH
