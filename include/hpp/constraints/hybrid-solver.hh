#ifndef HPP_CONSTRAINTS_HYBRID_SOLVER_HH
#define HPP_CONSTRAINTS_HYBRID_SOLVER_HH

#include <vector>

#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/config.hh>

#include <hpp/constraints/explicit-solver.hh>
#include <hpp/constraints/iterative-solver.hh>

namespace hpp {
  namespace constraints {
    class HPP_CONSTRAINTS_DLLAPI HybridSolver
      : public HierarchicalIterativeSolver
    {
      public:
        HybridSolver (const std::size_t& argSize, const std::size_t derSize)
          : HierarchicalIterativeSolver(), explicit_ (argSize, derSize),
          JeExpanded_ (derSize, derSize)
        {}

        ExplicitSolver& explicitSolver()
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

        template <typename LineSearchType>
        bool oneStep (vectorIn_t arg, LineSearchType& lineSearch) const
        {
          computeValue<true> (arg);
          updateJacobian (arg);
          computeDescentDirection ();
          lineSearch (*this, arg, dq_);
          explicit_.solve (arg);
          return isSatisfied(arg);
        }

      private:
        template <typename LineSearchType>
        Status impl_solve (vectorOut_t arg, LineSearchType ls) const;

        void updateJacobian (vectorIn_t arg) const;

        ExplicitSolver explicit_;
        mutable matrix_t Je_, JeExpanded_;
    }; // class HybridSolver
  } // namespace constraints
} // namespace hpp

#include <hpp/constraints/impl/hybrid-solver.hh>

#endif // HPP_CONSTRAINTS_HYBRID_SOLVER_HH
