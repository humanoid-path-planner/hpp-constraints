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

#ifndef HPP_CONSTRAINTS_EXPLICIT_CONSTRAINT_SET_HH
#define HPP_CONSTRAINTS_EXPLICIT_CONSTRAINT_SET_HH

#include <vector>

#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/config.hh>

#include <hpp/constraints/matrix-view.hh>
#include <hpp/constraints/differentiable-function-set.hh>

namespace hpp {
  namespace constraints {
    /// \addtogroup solvers
    /// \{

    /**
    Set of explicit constraints

    This class combines compatible explicit constraints as
    defined in the following paper published in Robotics Science and System
    2018: https://hal.archives-ouvertes.fr/hal-01804774/file/paper.pdf,
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

    To add explicit constraints, use methods ExplicitConstraintSet::add. If the
    constraint to add is not compatible with the previous one, this method
    returns -1.

    Method ExplicitConstraintSet::solve solves the explicit constraints.

    The combination of compatible explicit constraints is an explicit
    constraint. As such this class can be considered as an explicit constraint.

    We will therefore use the following notation
    \li \f$in\f$ for the set of indices of input variables,
    \li \f$out\f$ for the set of indices of output variables.
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

        /// Size of error vector
        /// \note that this size may differ from size of \link
        ///       ExplicitConstraintSet::outDers outDers output\endlink,
        size_type errorSize() const
        {
          return errorSize_;
        }

        /// Whether input vector satisfies the constraints of the solver
        /// \param arg input vector
	/// \param errorThreshold threshold to compare against the norms of
	///        the constraint errors. Default value is accessible via
	///        methods errorThreshold.
        bool isSatisfied (vectorIn_t arg, value_type errorThreshold = -1) const;

        /// Whether input vector satisfies the constraints of the solver
        /// \param arg input vector
        /// \retval error the constraint errors
	/// \param errorThreshold threshold to compare against the norms of
	///        the constraint errors. Default value is accessible via
	///        methods errorThreshold.
        bool isSatisfied (vectorIn_t arg, vectorOut_t error,
			  value_type errorThreshold = -1) const;

        /// Whether a constraint is satisfied for an input vector
        ///
        /// \param constraint, the constraint in the solver,
        /// \param arg the input vector
        /// \retval error the error of the constraint.
        /// \retval constraintFound whether the constraint belongs to the
        ///         solver,
        /// \return true if constraint belongs to solver and error is below
        ///         the threshold, false otherwise.
        bool isConstraintSatisfied (const ImplicitPtr_t& constraint,
                                    vectorIn_t arg, vectorOut_t error,
                                    bool& constraintFound) const;
        /// \}

        /// \name Construction of the problem
        /// \{

        /// Attempt to add an explicit constraint
        ///
        /// \param constraint explicit constraint
        /// \return the index of the function if the function was added,
        /// -1 otherwise.
        /// \note A function can be added iff it is compatible with the
        ///       previously added functions.
        size_type add (const ExplicitPtr_t& constraint);

        /// Check whether an explicit numerical constraint has been added
        /// \param numericalConstraint explicit numerical constraint
        /// \return true if the constraint is in the set.
        /// \note Comparison between constraints is performed by
        /// function names. This means that two constraints with the
        /// same function names are considered as equal.
        bool contains (const ExplicitPtr_t& numericalConstraint) const;

        /// Constructor
        ///
        /// \param space Lie group on which constraints are defined.
        ExplicitConstraintSet (const LiegroupSpacePtr_t& space)
          : configSpace_ (space)
          ,   inArgs_ (), notOutArgs_ ()
          ,   inDers_ (), notOutDers_ ()
          ,  outArgs_ (),  outDers_ ()
          , argFunction_ (Eigen::VectorXi::Constant(space->nq (), -1))
          , derFunction_ (Eigen::VectorXi::Constant(space->nv (), -1))
          , errorThreshold_ (Eigen::NumTraits<value_type>::epsilon())
          , errorSize_(0)
          // , Jg (nv, nv)
          , arg_ (space->nq ()), diff_(space->nv ()), diffSmall_()
        {
          notOutArgs_.addRow(0, space->nq ());
          notOutDers_.addCol(0, space->nv ());
        }

        /// \}

        /// \name Parameters
        /// \{

        /// Set error threshold
        void errorThreshold (const value_type& threshold)
        {
          errorThreshold_ = threshold;
        }
        /// Get error threshold
        value_type errorThreshold () const
        {
          return errorThreshold_;
        }
        /// Get error threshold
        value_type squaredErrorThreshold () const
        {
          return errorThreshold_*errorThreshold_;
        }

        /// \}

        /// \name Input and outputs
        /// \{

        /// Set \f$in\f$ of input configuration variables
        const RowBlockIndices& inArgs () const
        {
          return inArgs_;
        }

        /// Set of input velocity variables
        const ColBlockIndices& inDers () const
        {
          return inDers_;
        }

        /// Set \f$i\overline{n\cup ou}t\f$ of other configuration variables
        ///
        /// Configuration variables that are not output variables
        const RowBlockIndices& notOutArgs () const
        {
          return notOutArgs_;
        }

        /// Set \f$i\overline{n\cup ou}t\f$ of other velocity variables
        ///
        /// Velocity variables that are not output variables
        const ColBlockIndices& notOutDers () const
        {
          return notOutDers_;
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

        /// Set \f$out\f$ of output configuration variables
        /// \return the set of intervals corresponding the the configuration
        ///         variables that are ouputs of the combination of explicit
        ///         constraints.
        const RowBlockIndices& outArgs () const
        {
          return outArgs_;
        }

        /// Set of output velocity variables
        /// \return the set of intervals corresponding the the velocity
        ///         variables that are ouputs of the combination of explicit
        ///         constraints.
        const RowBlockIndices& outDers () const
        {
          return outDers_;
        }

        /// The Lie group on which constraints are defined
        LiegroupSpacePtr_t configSpace () const
        {
          return configSpace_;
        }

        /// The number of variables
        std::size_t nq () const
        {
          return configSpace_->nq ();
        }

        /// The number of derivative variables
        std::size_t nv () const
        {
          return configSpace_->nv ();
        }

        /// \}

        /// Return Jacobian matrix of output variable wrt not output variables
        ///
        /// \retval jacobian Jacobian matrix of the mapping from non output
        ///         variables to output variables. The columns of this matrix
        ///         corresponding to variables \f$in\f$ are filled with the
        ///         Jacobian of \f$f\f$:
        ///         \f{equation}
        ///         \frac{\partial f}{\partial \mathbf{q}_{in}}
        ///         (\mathbf{q}_{in}).
        ///         \f}
        ///         The columns corresponding to variables
        ///         \f$i\overline{n\cup ou}t\f$ are set to 0.
        inline MatrixBlockView jacobianNotOutToOut (matrix_t& jacobian) const
        {
          return MatrixBlockView(jacobian,
              outDers_.nbIndices() , outDers_.indices(),
              notOutDers_.nbIndices(), notOutDers_.indices());
        }

        /** Compute the Jacobian of the explicit constraint resolution

            \param q input configuration
            \param jacobian square Jacobian matrix of same size as velocity
                            i.e. given by \ref nv method.

            The result is the Jacobian of the explicit constraint set considered
            as a projector that maps to any \f$\mathbf{p}\in\mathcal{C}\f$,
            \f$\mathbf{q} = E(\mathbf{p})\f$ defined by
            \f{eqnarray}
            \mathbf{q}_{\bar{out}} &=& \mathbf{p}_{\bar{out}} \\
            \mathbf{q}_{out} &=& f (\mathbf{p}_{in}) + rhs
            \f}

            \warning it is assumed solve(q) has been called before.
        */
        void jacobian(matrixOut_t jacobian, vectorIn_t q) const;

        /// \name Right hand side accessors
        /// \{

        /// Compute right hand side of constraints using input configuration.
        ///
        /// \param p vector in \f$\mathcal{C}\f$.
        ///
        /// For each explicit constraint \f$E=(in,out,f)\f$, compute the right
        /// hand side as follows:
        /// \f{equation}
        /// rhs = f (\mathbf{q}),
        /// \f}
        /// where in general
        ///\f{equation}
        /// f(\mathbf{q}) = \mathbf{p}_{out} - f(\mathbf{q}_{in),
        /// \f}
        /// in such a way that all \f$\mathbf{q}\f$ satisfies the explicit
        /// constraint.
        /// \note For hpp::constraints::explicit_::RelativePose, the implicit
        ///       formulation does not derive from the explicit one. The
        ///       right hand side considered is the right hand side of the
        ///       implicit formulation.
        vector_t rightHandSideFromInput (vectorIn_t p);

        /// Compute right hand side of constraint using input configuration.
        ///
        /// \param constraint explicit constraint,
        /// \param p vector in \f$\mathcal{C}\f$.
        ///
        /// Let \f$E=(in,out,f)\f$ be the explicit constraint, compute the right
        /// hand side as follows:
        /// \f{equation}
        /// rhs = f (\mathbf{q}),
        /// \f}
        /// where in general
        ///\f{equation}
        /// f(\mathbf{q}) = \mathbf{p}_{out} - f(\mathbf{q}_{in),
        /// \f}
        /// in such a way that all \f$\mathbf{q}\f$ satisfies the explicit
        /// constraint.
        /// \note For hpp::constraints::explicit_::RelativePose, the implicit
        ///       formulation does not derive from the explicit one. The
        ///       right hand side considered is the right hand side of the
        ///       implicit formulation.
        bool rightHandSideFromInput (const ExplicitPtr_t& constraint,
                                     vectorIn_t p);

        /// Compute right hand side of constraint using input configuration.
        ///
        /// \param i index of the explicit constraint,
        /// \param p vector in \f$\mathcal{C}\f$.
        ///
        /// Let \f$E=(in,out,f)\f$ be the explicit constraint, compute the right
        /// hand side as follows:
        /// \f{equation}
        /// rhs = f (\mathbf{q}),
        /// \f}
        /// where in general
        ///\f{equation}
        /// f(\mathbf{q}) = \mathbf{p}_{out} - f(\mathbf{q}_{in),
        /// \f}
        /// in such a way that all \f$\mathbf{q}\f$ satisfies the explicit
        /// constraint.
        /// \note For hpp::constraints::explicit_::RelativePose, the implicit
        ///       formulation does not derive from the explicit one. The
        ///       right hand side considered is the right hand side of the
        ///       implicit formulation.
        void rightHandSideFromInput (const size_type& i, vectorIn_t p);

        /// Set the right hand sides of the explicit constraints.
        ///
        /// \param rhs the right hand side.
        ///
        /// The components of rhs are dispatched to the right hand sides of the
        /// explicit constraints in the order they are added.
        void rightHandSide (vectorIn_t rhs);

        /// Set the right hand side for a given explicit constraint
        ///
        /// \param constraint the explicit constraint,
        /// \param rhs right hand side.
        bool rightHandSide (const ExplicitPtr_t& constraint, vectorIn_t rhs);

	/// Get the right hand side for a given explicit constraint
        ///
        /// \param constraint the explicit constraint,
        /// \param rhs right hand side.
        bool getRightHandSide (const ExplicitPtr_t& constraint, vectorOut_t rhs) const ;

        /// Set the right hand side for a given explicit constraint
        ///
        /// \param i order of the explicit constraint,
        /// \param rhs right hand side.
        void rightHandSide (const size_type& i, vectorIn_t rhs);

        /// Get the right hand sides
        /// \return the right hand sides of the explicit constraints stacked
        ///         into a vector
        vector_t rightHandSide () const;

        /// Get size of the right hand side
        ///
        /// See documentation of classes Implicit and Explicit for details.
        size_type rightHandSideSize () const;

        /// \}

        std::ostream& print (std::ostream& os) const;

      private:
        typedef std::vector<bool> Computed_t;

        /// Compute output variables with respect to input variables
        /// \param i index of explicit constraint,
        /// \retval arg configuration of the system in which output variables
        ///             are set to their values.
        void solveExplicitConstraint(const std::size_t& i, vectorOut_t arg)
          const;
        /// Compute rows of Jacobian corresponding to output of function
        ///
        /// \param i index of the explicit constraint,
        /// \retval J Jacobian in which rows are computed
        ///
        /// Let
        ///   \li E = (f, in, out) be the explicit constraint of index i,
        ///   \li E.jacobian be the Jacobian of f,
        ///   \li E.in the input velocity variables of the constraints,
        ///   \li E.out the output velocity variables of the constraints,
        ///   \li Jin the matrix composed of E.in rows of J,
        ///   \li Jout the matrix composed of E.out rows of J,
        /// then,
        ///   Jout = E.jacobian * Jin
        void computeJacobian(const std::size_t& i, matrixOut_t J) const;
        void computeOrder(const std::size_t& iF, std::size_t& iOrder, Computed_t& computed);

        LiegroupSpacePtr_t configSpace_;

        struct Data {
          Data (const ExplicitPtr_t& constraint);
          ExplicitPtr_t constraint;
          RowBlockIndices equalityIndices;
          LiegroupElement rhs_implicit;
          // implicit formulation
          mutable LiegroupElement h_value;
          // explicit formulation
          mutable vector_t qin, qout;
          mutable LiegroupElement f_value, res_qout;
          // jacobian of f
          mutable matrix_t jacobian;
        }; // struct Data

        RowBlockIndices inArgs_, notOutArgs_;
        ColBlockIndices inDers_, notOutDers_;
        /// Output indices
        RowBlockIndices outArgs_, outDers_;

        Eigen::MatrixXi inOutDependencies_;

        std::vector<Data> data_;
        std::vector<std::size_t> computationOrder_;
        /// For each configuration variable i, argFunction_[i] is the index in
        /// data_ of the function that computes this configuration
        /// variable.
        /// -1 means that the configuration variable is not ouput of any
        /// function in data_.
        Eigen::VectorXi argFunction_, derFunction_;
        value_type errorThreshold_;
        size_type errorSize_;
        // mutable matrix_t Jg;
        mutable vector_t arg_, diff_, diffSmall_;

        /// Constructor for serialization
        ExplicitConstraintSet() 
          :  inArgs_ (), notOutArgs_ ()
          ,  inDers_ (), notOutDers_ ()
          , outArgs_ (),  outDers_ ()
          , errorThreshold_ (Eigen::NumTraits<value_type>::epsilon())
          , errorSize_(0)
        {}
        /// Initialization for serialization
        void init(const LiegroupSpacePtr_t& space)
        {
          configSpace_ = space;
          argFunction_ = Eigen::VectorXi::Constant(space->nq (), -1);
          derFunction_ = Eigen::VectorXi::Constant(space->nv (), -1);
          arg_.resize(space->nq ());
          diff_.resize(space->nv ());

          notOutArgs_.addRow(0, space->nq ());
          notOutDers_.addCol(0, space->nv ());
        }
        
        friend class solver::BySubstitution;
    }; // class ExplicitConstraintSet
    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_EXPLICIT_CONSTRAINT_SET_HH
