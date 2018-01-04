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

#ifndef HPP_CONSTRAINTS_HYBRID_SOLVER_HH
#define HPP_CONSTRAINTS_HYBRID_SOLVER_HH

#include <vector>

#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/config.hh>

#include <hpp/constraints/explicit-solver.hh>
#include <hpp/constraints/iterative-solver.hh>

namespace hpp {
  namespace constraints {
    /// \addtogroup solvers
    /// \{
    class HPP_CONSTRAINTS_DLLAPI HybridSolver
      : public HierarchicalIterativeSolver
    {
      public:
        HybridSolver (const std::size_t& argSize, const std::size_t derSize)
          : HierarchicalIterativeSolver(argSize, derSize), explicit_ (argSize, derSize),
          JeExpanded_ (derSize, derSize)
        {}

        ExplicitSolver& explicitSolver()
        {
          return explicit_;
        }

        const ExplicitSolver& explicitSolver () const
        {
          return explicit_;
        }

        /// Should be called whenever explicit solver is modified
        void explicitSolverHasChanged();

        template <typename LineSearchType>
        Status solve (vectorOut_t arg, LineSearchType ls = LineSearchType()) const
        {
          // TODO when there are only locked joint explicit constraints,
          // there is no need for this intricated loop.
          // if (explicit_.isConstant()) {
            // explicit_.solve(arg);
            // iterative_.solve(arg, ls);
          // } else {
          return impl_solve (arg, ls);
          // }
        }

        inline Status solve (vectorOut_t arg) const
        {
          return solve(arg, DefaultLineSearch());
        }

        bool isSatisfied (vectorIn_t arg) const
        {
          return 
            HierarchicalIterativeSolver::isSatisfied (arg)
            && explicit_.isSatisfied (arg);
        }

        bool isSatisfied (vectorIn_t arg, vectorOut_t error) const
        {
          assert (error.size() == dimension() + explicit_.outDers().nbIndices());
          bool iterative =
            HierarchicalIterativeSolver::isSatisfied (arg);
          residualError(error.head(dimension()));
          bool _explicit =
            explicit_.isSatisfied (arg, error.tail(explicit_.outDers().nbIndices()));
          return iterative && _explicit;
        }

        /// Project the point arg + darg onto the null space of the jacobian
        /// at arg.
        void projectOnKernel (vectorIn_t arg, vectorIn_t darg, vectorOut_t result) const;

        template <typename LineSearchType>
        bool oneStep (vectorOut_t arg, LineSearchType& lineSearch) const
        {
          computeValue<true> (arg);
          updateJacobian (arg);
          computeDescentDirection ();
          lineSearch (*this, arg, dq_);
          explicit_.solve (arg);
          return HierarchicalIterativeSolver::isSatisfied(arg);
        }

        /// Computes the jacobian of the explicit functions and
        /// updates the jacobian of the problem using the chain rule.
        void updateJacobian (vectorIn_t arg) const;

        /// Set error threshold
        void errorThreshold (const value_type& threshold)
        {
          HierarchicalIterativeSolver::errorThreshold(threshold);
          explicit_.errorThreshold(threshold);
        }
        /// Get error threshold
        value_type errorThreshold () const
        {
          return HierarchicalIterativeSolver::errorThreshold();
        }

        /// \name Right hand side accessors
        /// \{

        /// Compute a right hand side using the input arg.
        vector_t rightHandSideFromInput (vectorIn_t arg)
        {
          const size_type top = parent_t::rightHandSideSize();
          const size_type bot = explicit_.rightHandSideSize();
          vector_t rhs (top + bot);
          rhs.head(top) = parent_t::rightHandSideFromInput (arg);
          rhs.tail(bot) = explicit_.rightHandSideFromInput (arg);
          return rhs;
        }

        /// Set the right hand side for a given constraint.
        /// \param fImplicit implicit formulation of the constraint. Can be NULL
        /// \param fExplicit explicit formulation of the constraint. Can be NULL
        /// \param arg a vector of size argSize_
        /// \warning At least one of fImplicit and fExplicit must be non-NULL.
        bool rightHandSideFromInput (
            const DifferentiableFunctionPtr_t& fImplicit,
            const DifferentiableFunctionPtr_t& fExplicit,
            vectorIn_t arg)
        {
          assert (fImplicit || fExplicit);
          if (fExplicit && explicit_.rightHandSideFromInput (fExplicit, arg))
            return true;
          if (fImplicit && parent_t::rightHandSideFromInput (fImplicit, arg))
            return true;
          return false;
        }

        /// Set the right hand side for a given constraint.
        /// \param fImplicit implicit formulation of the constraint. Can be NULL
        /// \param fExplicit explicit formulation of the constraint. Can be NULL
        /// \param rhs the desired right hand side
        /// \warning At least one of fImplicit and fExplicit must be non-NULL.
        bool rightHandSide (
            const DifferentiableFunctionPtr_t& fImplicit,
            const DifferentiableFunctionPtr_t& fExplicit,
            vectorIn_t rhs)
        {
          assert (fImplicit || fExplicit);
          if (fExplicit && explicit_.rightHandSide (fExplicit, rhs))
            return true;
          if (fImplicit && parent_t::rightHandSide (fImplicit, rhs))
            return true;
          return false;
        }

        /// Set the level set parameter.
        /// \param rhs the level set parameter.
        void rightHandSide (vectorIn_t rhs)
        {
          const size_type top = parent_t::rightHandSideSize();
          const size_type bot = explicit_.rightHandSideSize();
          parent_t::rightHandSide (rhs.head(top));
          explicit_.rightHandSideFromInput (rhs.head(bot));
        }

        /// Get the level set parameter.
        /// \return the parameter.
        vector_t rightHandSide () const
        {
          const size_type top = parent_t::rightHandSideSize();
          const size_type bot = explicit_.rightHandSideSize();
          vector_t rhs (top + bot);
          rhs.head(top) = parent_t::rightHandSide ();
          rhs.tail(bot) = explicit_.rightHandSide ();
          return rhs;
        }

        /// Get size of the level set parameter.
        size_type rightHandSideSize () const
        {
          const size_type top = parent_t::rightHandSideSize();
          const size_type bot = explicit_.rightHandSideSize();
          return top + bot;
        }

        /// \}

        virtual std::ostream& print (std::ostream& os) const;

      protected:
        void integrate(vectorIn_t from, vectorIn_t velocity, vectorOut_t result) const
        {
          HierarchicalIterativeSolver::integrate(from, velocity, result);
          explicit_.solve (result);
        }

        void computeActiveRowsOfJ (std::size_t iStack);

      private:
        typedef HierarchicalIterativeSolver parent_t;

        template <typename LineSearchType>
        Status impl_solve (vectorOut_t arg, LineSearchType ls) const;

        ExplicitSolver explicit_;
        mutable matrix_t Je_, JeExpanded_;

        friend struct lineSearch::Backtracking;
    }; // class HybridSolver
    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_HYBRID_SOLVER_HH
