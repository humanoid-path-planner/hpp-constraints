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

#include <boost/function.hpp>

#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/config.hh>

#include <hpp/constraints/matrix-view.hh>
#include <hpp/constraints/differentiable-function-stack.hh>

namespace hpp {
  namespace constraints {
    /// \addtogroup solvers
    /// \{

    /// Solve system of explicit functions
    ///
    /// The solver works on a given set of variables \f$ X = (x_i) \f$.
    /// It contains a set of functions \f$ f_j \f$ that takes as input a subset of \f$ X \f$ and
    /// outputs values corresponding to another subset of \f$ X \f$.
    /// The equation are then of the following form:
    /// \f$ X_{output} = f (X_{input}) + rhs \f$ where the addition is the
    /// addition of the \ref outputSpace() of \f$ f \f$.
    /// 
    /// There can be no cycles in the dependencies. Moreover, a variable cannot
    /// be the output of two different functions.
    /// For instance, \f$ (x_0, x_2) = f_0( x_1, x_3 ) \f$ and
    /// \f$ (x_3) = f_1( x_4 ) \f$ is a valid input. It would not be possible to
    /// add \f$ (x_0) = f_2( x_2 ) \f$ because it would introduce a cycle,
    /// or \f$ (x_3) = f_3( x_1 ) \f$ because \f$ x_3 \f$ would be computed by
    /// two different functions.
    ///
    /// The resolution consists in modyfing the output values of each function,
    /// while respecting the dependendy order. Considering \f$ f_0, f_1 \f$
    /// above, \f$ f_1 \f$ must be computed first.
    class HPP_CONSTRAINTS_DLLAPI ExplicitSolver
    {
      public:
        typedef Eigen::RowBlockIndices RowBlockIndices;
        typedef Eigen::ColBlockIndices ColBlockIndices;
        typedef Eigen::MatrixBlockView<matrix_t, Eigen::Dynamic, Eigen::Dynamic, false, false> MatrixBlockView;

        /// \name Resolution
        /// \{

        bool solve (vectorOut_t arg) const;

        bool isSatisfied (vectorIn_t arg) const;

        bool isSatisfied (vectorIn_t arg, vectorOut_t error) const;
        
        /// \}

        /// \name Construction of the problem
        /// \{

        /// Returns true if the function was added, false otherwise
        /// A function can be added iif its outputs do not interfere with the
        /// output of another function.
        bool add (const DifferentiableFunctionPtr_t& f,
            const RowBlockIndices& inArg,
            const RowBlockIndices& outArg,
            const ColBlockIndices& inDer,
            const RowBlockIndices& outDer);

        /// Returns true if the function was added, false otherwise
        /// A function can be added iif its outputs do not interfere with the
        /// output of another function.
        bool add (const DifferentiableFunctionPtr_t& f,
            const RowBlockIndices& inArg,
            const RowBlockIndices& outArg,
            const ColBlockIndices& inDer,
            const RowBlockIndices& outDer,
            const ComparisonTypes_t& comp);

        /// \warning the two functions must have the same input and output
        /// indices.
        bool replace (const DifferentiableFunctionPtr_t& oldf,
                      const DifferentiableFunctionPtr_t& newd);

        ExplicitSolver (const std::size_t& argSize, const std::size_t derSize)
          : argSize_ (argSize), derSize_ (derSize)
          ,   inArgs_ (), freeArgs_ ()
          ,   inDers_ (), freeDers_ ()
          ,  outArgs_ (),  outDers_ ()
          , argFunction_ (Eigen::VectorXi::Constant(argSize, -1))
          , derFunction_ (Eigen::VectorXi::Constant(derSize, -1))
          , squaredErrorThreshold_ (Eigen::NumTraits<value_type>::epsilon())
          // , Jg (derSize, derSize)
          , arg_ (argSize), diff_(derSize), diffSmall_()
        {
          freeArgs_.addRow(0, argSize);
          freeDers_.addCol(0, derSize);
        }

        /// \}

        /// \name Parameters
        /// \{

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

        /// \}

        /// \name Input and outputs
        /// \{

        /// The set of variable indices which affects the output.
        /// This is a subset of \ref freeArgs
        const RowBlockIndices& inArgs () const
        {
          return inArgs_;
        }

        /// The set of derivative variable indices which affects the output.
        /// This is a subset of \ref freeDers
        const ColBlockIndices& inDers () const
        {
          return inDers_;
        }

        /// The set of variable indices which are not affected by the
        /// resolution.
        const RowBlockIndices& freeArgs () const
        {
          return freeArgs_;
        }

        /// The set of derivative variable indices which are not affected by the
        /// resolution.
        const ColBlockIndices& freeDers () const
        {
          return freeDers_;
        }

        /// Same as \ref inArgs
        ColBlockIndices activeParameters () const;

        /// Same as \ref inDers
        const ColBlockIndices& activeDerivativeParameters () const;

        /// Returns a matrix of integer whose:
        /// - rows correspond to functions
        /// - cols correspond to DoF
        /// - values correspond to the dependency degree of a function wrt to
        ///   a DoF
        const Eigen::MatrixXi& inOutDependencies () const
        {
          return inOutDependencies_;
        }

        /// Same as \ref inOutDependencies except that cols correpond to DoFs.
        Eigen::MatrixXi inOutDofDependencies () const;

        const Eigen::VectorXi& derFunction () const
        {
          return derFunction_;
        }

        /// The set of variable indices which are computed.
        const RowBlockIndices& outArgs () const
        {
          return outArgs_;
        }

        /// The set of derivative variable indices which are computed.
        const RowBlockIndices& outDers () const
        {
          return outDers_;
        }

        /// The number of variables
        const std::size_t& argSize () const
        {
          return argSize_;
        }

        /// The number of derivative variables
        const std::size_t& derSize () const
        {
          return derSize_;
        }

        /// \}

        inline MatrixBlockView viewJacobian(matrix_t& jacobian) const
        {
          return MatrixBlockView(jacobian,
              outDers_.nbIndices() , outDers_.indices(),
              freeDers_.nbIndices(), freeDers_.indices());
        }

        // /// \param jacobian must be of dimensions (derSize - freeDers().nbIndices(), freeDers().nbIndices())
        /// \param jacobian must be of dimensions (derSize, derSize) but only a subsegment will be used.
        /// \warning it is assumed solve(arg) has been called before.
        void jacobian(matrixOut_t jacobian, vectorIn_t arg) const;

        /// \name Right hand side accessors
        /// \{

        /// Compute a right hand side using the input arg.
        vector_t rightHandSideFromInput (vectorIn_t arg);

        bool rightHandSideFromInput (const DifferentiableFunctionPtr_t& f, vectorIn_t arg);

        /// Set the level set parameter.
        /// \param rhs the level set parameter.
        void rightHandSide (vectorIn_t rhs);

        bool rightHandSide (const DifferentiableFunctionPtr_t& f, vectorIn_t rhs);

        /// Get the level set parameter.
        /// \return the parameter.
        vector_t rightHandSide () const;

        /// Get size of the level set parameter.
        size_type rightHandSideSize () const;

        /// \}

        std::ostream& print (std::ostream& os) const;

      private:
        typedef std::vector<bool> Computed_t;

        void computeFunction(const std::size_t& i, vectorOut_t arg) const;
        void computeJacobian(const std::size_t& i, matrixOut_t J) const;
        void computeOrder(const std::size_t& iF, std::size_t& iOrder, Computed_t& computed);

        const std::size_t argSize_, derSize_;

        struct Function {
          Function (DifferentiableFunctionPtr_t _f, RowBlockIndices ia,
                    RowBlockIndices oa, ColBlockIndices id, RowBlockIndices od,
                    const ComparisonTypes_t& comp);
          void setG (const DifferentiableFunctionPtr_t& _g, const DifferentiableFunctionPtr_t& _ginv);
          DifferentiableFunctionPtr_t f;
          DifferentiableFunctionPtr_t g, ginv;
          RowBlockIndices inArg, outArg;
          ColBlockIndices inDer;
          RowBlockIndices outDer;
          ComparisonTypes_t comparison;
          RowBlockIndices equalityIndices;
          vector_t rightHandSide;

          mutable vector_t qin, qout;
          mutable LiegroupElement value, expected;
          mutable matrix_t jacobian, jGinv;
        }; // struct Function

        RowBlockIndices inArgs_, freeArgs_;
        ColBlockIndices inDers_, freeDers_;
        RowBlockIndices outArgs_, outDers_;

        Eigen::MatrixXi inOutDependencies_;

        std::vector<Function> functions_;
        std::vector<std::size_t> computationOrder_;
        /// For dof i, dofFunction_[i] is the index of the function that computes it.
        /// -1 means it is the output of no function.
        Eigen::VectorXi argFunction_, derFunction_;
        value_type squaredErrorThreshold_;
        // mutable matrix_t Jg;
        mutable vector_t arg_, diff_, diffSmall_;
    }; // class ExplicitSolver
    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_EXPLICIT_SOLVER_HH
