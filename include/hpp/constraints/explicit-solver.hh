
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

#ifndef HPP_CONSTRAINTS_EXPLICIT_SOLVER_HH
#define HPP_CONSTRAINTS_EXPLICIT_SOLVER_HH

#include <deque>

#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/config.hh>

#include <hpp/constraints/matrix-view.hh>
#include <hpp/constraints/differentiable-function-stack.hh>

namespace hpp {
  namespace constraints {
    /*
    struct HPP_CONSTRAINTS_DLLAPI ExplicitFunction :
      public DifferentiableFunction
    {
      public:
        ExplicitFunction (const DifferentiableFunctionPtr_t& function,
	 const SizeIntervals_t& outputConf,
	 const SizeIntervals_t& outputVelocity)
          : function_ (function)
        {}

      protected:
        virtual void impl_compute (vectorOut_t result,
                                   vectorIn_t argument) const
        {
          implicit_->value(value, argument);
          input_.view(argument).writeTo(result);
        }

        virtual void impl_jacobian (matrixOut_t jacobian,
                                    vectorIn_t arg) const
        {
        }

      private:
        typedef Eigen::MatrixBlockIndexes<false, true> RowBlockIndexes;

        DifferentiableFunctionPtr_t implicit_;
        RowBlockIndexes input_ ;
        RowBlockIndexes output_;
        mutable vector_t value;
        mutable matrix_t jacobian;
    }; // class ExplicitSolver
    */

    class HPP_CONSTRAINTS_DLLAPI ExplicitSolver
    {
      public:
        typedef Eigen::RowBlockIndexes RowBlockIndexes;

        bool solve (vectorOut_t arg) const;

        /// Returns true if the function was added, false otherwise
        /// A function can be added iif its outputs do not interfere with the
        /// output of another function.
        bool add (const DifferentiableFunctionPtr_t& f,
            const RowBlockIndexes& input,
            const RowBlockIndexes& output);

        ExplicitSolver (const std::size_t& argSize, const std::size_t derSize)
          : argSize_ (argSize), derSize_ (derSize)
          , freeDofs_ ()
          , dofFunction_ (Eigen::VectorXi::Constant(argSize, -1))
        {
          freeDofs_.addRow(0, argSize);
        }

        const RowBlockIndexes& freeDofs () const;

      private:
        typedef std::vector<bool> Computed_t;

        void computeFunction(const std::size_t& i, vectorOut_t arg, Computed_t& computed) const;

        const std::size_t argSize_, derSize_;

        struct Function {
          Function (DifferentiableFunctionPtr_t _f, RowBlockIndexes ia, RowBlockIndexes oa)
            : f (_f), inArg (ia), outArg (oa)
          {}
          DifferentiableFunctionPtr_t f;
          RowBlockIndexes inArg, outArg;
        };

        RowBlockIndexes freeDofs_;

        std::deque<Function> functions_;
        /// For dof i, dofFunction_[i] is the index of the function that computes it.
        /// -1 means it is the output of no function.
        Eigen::VectorXi dofFunction_;
        // /// If dof[i] depends on dof[j], then dof[j] must be computed before dof[i]
        // std::vector<int> dofDependency_;
        // /// The dof that 
        // std::vector<int> dofLeaves_;
    }; // class ExplicitSolver
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_EXPLICIT_SOLVER_HH
