
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

#include <vector>

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
        typedef Eigen::ColBlockIndexes ColBlockIndexes;
        // typedef Eigen::MatrixBlockIndexes<false, false> MatrixBlockIndexes;
        typedef Eigen::MatrixBlockView<matrix_t, Eigen::Dynamic, Eigen::Dynamic, false, false> MatrixBlockView;

        bool solve (vectorOut_t arg) const;

        /// Returns true if the function was added, false otherwise
        /// A function can be added iif its outputs do not interfere with the
        /// output of another function.
        bool add (const DifferentiableFunctionPtr_t& f,
            const RowBlockIndexes& inArg,
            const RowBlockIndexes& outArg,
            const ColBlockIndexes& inDer,
            const RowBlockIndexes& outDer);

        ExplicitSolver (const std::size_t& argSize, const std::size_t derSize)
          : argSize_ (argSize), derSize_ (derSize)
          ,  inArgs_ (),  inDers_ ()
          , outArgs_ (), outDers_ ()
          , argFunction_ (Eigen::VectorXi::Constant(argSize, -1))
          , derFunction_ (Eigen::VectorXi::Constant(derSize, -1))
          // , Jg (derSize, derSize)
        {
          inArgs_.addRow(0, argSize);
          inDers_.addCol(0, derSize);
        }

        const RowBlockIndexes& inArgs () const
        {
          return inArgs_;
        }

        const ColBlockIndexes& inDers () const
        {
          return inDers_;
        }

        const RowBlockIndexes& outArgs () const
        {
          return outArgs_;
        }

        const RowBlockIndexes& outDers () const
        {
          return outDers_;
        }

        const std::size_t& argSize () const
        {
          return argSize_;
        }

        const std::size_t& derSize () const
        {
          return derSize_;
        }

        inline MatrixBlockView viewJacobian(matrix_t& jacobian) const
        {
          return MatrixBlockView(jacobian,
              outDers_.nbIndexes(), outDers_.indexes(),
              inDers_.nbIndexes(), inDers_.indexes());
        }

        // /// \param jacobian must be of dimensions (derSize - freeDers().nbIndexes(), freeDers().nbIndexes())
        /// \param jacobian must be of dimensions (derSize, derSize) but only a subsegment will be used.
        /// \warning it is assumed solve(arg) has been called before.
        void jacobian(matrixOut_t jacobian, vectorIn_t arg) const;

      private:
        typedef std::vector<bool> Computed_t;

        void computeFunction(const std::size_t& i, vectorOut_t arg, Computed_t& computed) const;
        void computeJacobian(const std::size_t& i, matrixOut_t J) const;
        void computeOrder(const std::size_t& iF, std::size_t& iOrder, Computed_t& computed);

        const std::size_t argSize_, derSize_;

        struct Function {
          Function (DifferentiableFunctionPtr_t _f, RowBlockIndexes ia, RowBlockIndexes oa, ColBlockIndexes id, RowBlockIndexes od)
            : f (_f), inArg (ia), outArg (oa), inDer (id), outDer (od)
          {}
          DifferentiableFunctionPtr_t f;
          RowBlockIndexes inArg, outArg;
          ColBlockIndexes inDer;
          RowBlockIndexes outDer;

          mutable vector_t value;
          mutable matrix_t jacobian;
        };

        RowBlockIndexes inArgs_;
        ColBlockIndexes inDers_;
        RowBlockIndexes outArgs_, outDers_;

        std::vector<Function> functions_;
        std::vector<std::size_t> computationOrder_;
        /// For dof i, dofFunction_[i] is the index of the function that computes it.
        /// -1 means it is the output of no function.
        Eigen::VectorXi argFunction_, derFunction_;
        // mutable matrix_t Jg;
    }; // class ExplicitSolver
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_EXPLICIT_SOLVER_HH