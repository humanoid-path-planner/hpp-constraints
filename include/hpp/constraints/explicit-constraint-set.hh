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

#ifndef HPP_CONSTRAINTS_EXPLICIT_CONSTRAINT_SET_HH
#define HPP_CONSTRAINTS_EXPLICIT_CONSTRAINT_SET_HH

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

    /**
    Set of explicit constraints

    This class combines compatible explicit constraints as
    defined in the following paper published in Robotics Science and System
    2018: https://hal.archives-ouvertes.fr/hal-01804774/file/paper.pdf\endlink,
    Section II-B Definition 4.

    An explicit constraint \f$E=(in,out,f)\f$ on a robot
    configuration space \f$\mathcal{C}\f$ is defined by
    \li a subset of input indices
        \f$in\subset\{1,\cdots, \dim\mathcal{C}\}\f$,
    \li a subset of output indices
        \f$out\subset\{1,\cdots, \dim\mathcal{C}\}\f$,
    \li a smooth mapping \f$f\f$ from \f$\mathbf{R}^{|in|}\f$ to \f$\mathbf{R}^{|out|}\f$,
    satisfying the following properties:
    \li \f$in\cap out = \emptyset\f$,
    \li for any \f$\mathbf{p}\in\mathcal{C}\f$,
        \f$\mathbf{q} = E(\mathbf{p})\f$ is defined by
      \f{eqnarray}
        &\mathbf{q}_{\bar{out}} = \mathbf{p}_{\bar{out}}\\
        &\mathbf{q}_{out} = f (\mathbf{p}_{in}).
      \f}

    \note Right hand side.

    For manipulation planning, it is useful to handle a parameterizable
    right hand side \f$rhs\f$. The expression above thus becomes

    \f{equation}
    \mathbf{q}_{out} = f (\mathbf{p}_{in}) + rhs.
    \f}

    The right hand side may be set using the various methods
    ExplicitConstraintSet::rightHandSide and
    ExplicitConstraintSet::rightHandSideFromInput.

    \note For some applications like manipulation planning, an
    invertible function \f$ g \f$ (of known inverse \f$ g^{-1} \f$)
    can be specified for each explicit constraint \f$E\f$. The above expression
    then becomes:
    \f{equation}
    g(\mathbf{q}_{out}) = f(\mathbf{p}_{in}) + rhs
    \f}

    To add explicit constraints, use methods ExplicitConstraintSet::add. Note
    that explicit constraints should be compatible.

    Method ExplicitConstraintSet::solve solves the explicit constraints.
    **/
    class HPP_CONSTRAINTS_DLLAPI ExplicitConstraintSet
    {
      public:
        typedef Eigen::RowBlockIndices RowBlockIndices;
        typedef Eigen::ColBlockIndices ColBlockIndices;
        typedef Eigen::MatrixBlockView<matrix_t, Eigen::Dynamic, Eigen::Dynamic, false, false> MatrixBlockView;

        /// \name Resolution
        /// \{

        /// Solve the explicit constraints
        /// \param arg input configuration,
        /// \retval arg output configuration satisfying the explicit
        ///         constraints:
        ///         \f$\mathbf{q}_{out} = g^{-1}
        ///         \left(f(\mathbf{q}_{in}) + rhs\right)\f$
        bool solve (vectorOut_t arg) const;

        bool isSatisfied (vectorIn_t arg) const;

        bool isSatisfied (vectorIn_t arg, vectorOut_t error) const;

        /// \}

        /// \name Construction of the problem
        /// \{

        /// Attempt to add an explicit constraint
        ///
        /// \param f differentiable function,
        /// \param inArg subset of input indices,
        /// \param outArg subset of output indices,
        /// \param inDer subset of input indices for the velocity,
        /// \param outDer subset of output indices for the velocity,
        /// \return the index of the function if the function was added,
        /// -1 otherwise.
        /// \note A function can be added iff it is compatible with the
        ///       previously added functions.
        size_type add (const DifferentiableFunctionPtr_t& f,
            const RowBlockIndices& inArg,
            const RowBlockIndices& outArg,
            const ColBlockIndices& inDer,
            const RowBlockIndices& outDer);

        /// Attempt to add an explicit constraint
        ///
        /// \param f differentiable function,
        /// \param inArg subset of input indices,
        /// \param outArg subset of output indices,
        /// \param inDer subset of input indices for the velocity,
        /// \param outDer subset of output indices for the velocity,
        /// \param comp Comparison type,
        /// \return the index of the function if the function was added,
        /// -1 otherwise.
        /// \note A function can be added iff it is compatible with the
        ///       previously added functions.
        size_type add (const DifferentiableFunctionPtr_t& f,
            const RowBlockIndices& inArg,
            const RowBlockIndices& outArg,
            const ColBlockIndices& inDer,
            const RowBlockIndices& outDer,
            const ComparisonTypes_t& comp);

        /// Set \f$g\f$  and \f$g^{-1}\f$ functions
        bool setG (const DifferentiableFunctionPtr_t& f,
                   const DifferentiableFunctionPtr_t& g,
                   const DifferentiableFunctionPtr_t& ginv);

        /// \warning the two functions must have the same input and output
        /// indices.
        bool replace (const DifferentiableFunctionPtr_t& oldf,
                      const DifferentiableFunctionPtr_t& newd);

        /// Constructor
        ///
        /// \param argSize dimension of vector space in which the robot
        ///                configuration space is immersed.
        /// \param derSize dimension of tangent space to configuration space.
        ExplicitConstraintSet (const std::size_t& argSize, const std::size_t derSize)
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

        /// Compute right hand side of explicit constraints using input configuration.
        ///
        /// \param p vector in \f$\mathcal{C}\f$.
        ///
        /// For each explicit constraint \f$E=(in,out,f)\f$, compute the right
        /// hand side as follows:
        /// \f{equation}
        /// rhs = g (\mathbf{p}_{out}) - f(\mathbf{q}_{in})
        /// \f}
        /// in such a way that all \f$\mathbf{q}\f$ satisfies all the explicit
        /// constraints.
        vector_t rightHandSideFromInput (vectorIn_t p);

        /// Compute right hand side of explicit constraint using input configuration.
        ///
        /// \param f differentiable function associated to the explicit constraints
        /// \param p vector in \f$\mathcal{C}\f$.
        ///
        /// Let \f$E=(in,out,f)\f$ be the explicit constraint, compute the right
        /// hand side as follows:
        /// \f{equation}
        /// rhs = g (\mathbf{p}_{out}) - f(\mathbf{q}_{in})
        /// \f}
        /// in such a way that all \f$\mathbf{q}\f$ satisfies the explicit
        /// constraint.
        bool rightHandSideFromInput (const DifferentiableFunctionPtr_t& f, vectorIn_t p);

        /// Compute right hand side of explicit constraint using input configuration.
        ///
        /// \param fidx order of the explicit constraint,
        /// \param p vector in \f$\mathcal{C}\f$.
        ///
        /// Let \f$E=(in,out,f)\f$ be the explicit constraint, compute the right
        /// hand side as follows:
        /// \f{equation}
        /// rhs = g (\mathbf{p}_{out}) - f(\mathbf{q}_{in})
        /// \f}
        /// in such a way that all \f$\mathbf{q}\f$ satisfies the explicit
        /// constraint.
        void rightHandSideFromInput (const size_type& fidx, vectorIn_t p);

        /// Set the right hand sides of the explicit constraints.
        ///
        /// \param rhs the right hand side.
        ///
        /// The components of rhs are dispatched to the right hand sides of the
        /// explicit constraints in the order they are added.
        void rightHandSide (vectorIn_t rhs);

        /// Set the right hand side for a given explicit constraint
        ///
        /// \param f the differentiable function of the explicit constraint,
        /// \param rhs right hand side.
        bool rightHandSide (const DifferentiableFunctionPtr_t& f, vectorIn_t rhs);

        /// Set the right hand side for a given explicit constraint
        ///
        /// \param fidx order of the explicit constraint,
        /// \param rhs right hand side.
        void rightHandSide (const size_type& fidx, vectorIn_t rhs);

        /// Get the right hand sides
        /// \return the right hand sides of the explicit constraints stacked
        ///         into a vector
        vector_t rightHandSide () const;

        /// Get size of the right hand side
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
    }; // class ExplicitConstraintSet
    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_EXPLICIT_CONSTRAINT_SET_HH
