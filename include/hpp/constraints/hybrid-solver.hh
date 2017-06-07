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
          assert (error.size() == dimension() + explicit_.outDers().nbIndexes());
          bool iterative =
            HierarchicalIterativeSolver::isSatisfied (arg);
          residualError(error.head(dimension()));
          bool _explicit =
            explicit_.isSatisfied (arg, error.tail(explicit_.outDers().nbIndexes()));
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

      private:
        template <typename LineSearchType>
        Status impl_solve (vectorOut_t arg, LineSearchType ls) const;

        ExplicitSolver explicit_;
        mutable matrix_t Je_, JeExpanded_;
    }; // class HybridSolver
    /// \}
  } // namespace constraints
} // namespace hpp

#include <hpp/constraints/impl/hybrid-solver.hh>

#endif // HPP_CONSTRAINTS_HYBRID_SOLVER_HH
