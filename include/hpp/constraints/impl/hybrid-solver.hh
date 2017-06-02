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

#ifndef HPP_CONSTRAINTS_IMPL_HYBRID_SOLVER_HH
#define HPP_CONSTRAINTS_IMPL_HYBRID_SOLVER_HH

namespace hpp {
  namespace constraints {
    template <typename LineSearchType>
    inline HybridSolver::Status HybridSolver::impl_solve (
        vectorOut_t arg,
        LineSearchType lineSearch) const
    {
      hppDout (info, "before projection: " << arg.transpose ());
      assert (!arg.hasNaN());

      explicit_.solve(arg);

      size_type errorDecreased = 3, iter = 0;
      value_type previousSquaredNorm =
	std::numeric_limits<value_type>::infinity();

      // Fill value and Jacobian
      computeValue<true> (arg);
      computeError();

      hppDout (info, "squareNorm = " << squaredNorm_);

      while (squaredNorm_ > squaredErrorThreshold_ && errorDecreased &&
	     iter < maxIterations_) {

        // Update the jacobian using the jacobian of the explicit system.
        updateJacobian(arg);

        computeDescentDirection ();
        lineSearch (*this, arg, dq_);
        explicit_.solve(arg);

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

    inline void HybridSolver::projectOnKernel (vectorIn_t arg, vectorIn_t darg, vectorOut_t result) const
    {
      computeValue<true> (arg);
      updateJacobian(arg);
      getReducedJacobian (reducedJ_);

      svd_.compute (reducedJ_);

      reduction_.view(darg).writeTo(dqSmall_);

      vector_t tmp = getV1(svd_).adjoint() * dqSmall_;
      dqSmall_.noalias() -= getV1(svd_) * tmp;

      result = reduction_.view(dqSmall_);
    }
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_IMPL_HYBRID_SOLVER_HH
