// Copyright (c) 2015 - 2018 CNRS
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr), Florent Lamiraux
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


#ifndef HPP_CONSTRAINTS_IMPLICIT_HH
# define HPP_CONSTRAINTS_IMPLICIT_HH

# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/config.hh>
# include <hpp/constraints/comparison-types.hh>
# include <hpp/constraints/matrix-view.hh>

# include <hpp/util/serialization-fwd.hh>

# include <hpp/pinocchio/liegroup-element.hh>

namespace hpp {
  namespace constraints {
    /**
        \addtogroup constraints
        \{
    */
    /**
        This class represents a parameterizable numerical constraint that
        compares the output of a function \f$h\f$ to a right hand side Lie
	group element.

        \par Definition

        \li The function \f$h\f$ takes input in a configuration space
        \f$\mathcal{C}\f$ and output in a Lie group \f$\mathbf{L}\f$,
        \li the dimensions of \f$\mathbf{L}\f$ and of its tangent space are
        respectively \f$(n_q,n_v)\f$.
        \li The comparison is represented by a vector \f$\mathbf{c}\f$ of
        dimension \f$n_v\f$ with values in
        enum hpp::contraints::ComparisonType = {
        \f$\mathbf{Equality}\f$, \f$\mathbf{EqualToZero}\f$,
        \f$\mathbf{Inferior}\f$, \f$\mathbf{Superior}\f$ }.
        \li The right hand side is Lie group element of dimension \f$n_q\f$.

        \par Error

        A configuration \f$\mathbf{q}\f$ is said to satisfy the constraint for
        a given right hand side if and only if the error \f$e\f$ as computed
        below is
        smaller in norm than a given threshold.

        Let
        \f[
        \Delta = h (\mathbf{q}) - rhs \in \mathbf{R}^{n_v},
        \f]
        for each component \f$i\in\{0,\cdots,n_v-1\}\f$,
        \li if \f$c_i\f$ is \f$\mathbf{Inferior}\f$,
        \f$e_i = \max (0,\Delta_i)\f$,
        \li if \f$c_i\f$ is \f$\mathbf{Superior}\f$,
        \f$e_i = \min (0,\Delta_i)\f$,
        \li if \f$c_i\f$ is \f$\mathbf{Equality}\f$,
        \f$e_i = \Delta_i\f$,
        \li if \f$c_i\f$ is \f$\mathbf{EqualToZero}\f$, \f$e_i = \Delta_i\f$.

	\par Mask

	A mask is a vector of Boolean values of size \f$n_v\f$. Values set to
	false means that the corresponding component of the error defined above
	is not taken into account to determine whether the constraint is
	satisfied. The active rows of the constraint may be accessed via
	method activeRows.

        \par Parameterizable right hand side

        Lines with \textbf{Equality} comparator in the above definition of the
        error need a parameter, while lines with comparators \textbf{Inferior},
        \textbf{Superior}, or \textbf{EqualToZero} do not.
        As a consequence, the right hand side of the constraint is defined by
        a vector \f$\lambda\f$ of parameters of size the number of
        \textbf{Equality} occurences in vector \f$\mathbf{c}\f$. The right hand
        side is then defined as in the following example:
        \f[
        rhs = \exp\left(\begin{array}{c}\lambda_1 \\ 0 \\ 0 \\ \lambda_2 \\
        \vdots \end{array}\right)
         \ \ \ \ \mathbf{c} = \left(\begin{array}{c}\mathbf{Equality} \\
        \mathbf{EqualToZero} \\ \mathbf{Inferior} \\ \mathbf{Equality} \\
        \vdots \end{array}\right)
        \f]
        To retrieve the size of vector \f$\lambda\f$, call method
        Implicit::parameterSize (). To set and get the right hand side value,
        use method Implicit::rightHandSide.

        \par Time varying right hand side

        The right hand side of the constraint may depend on time, for instance
        if the constraint is associated to a trajectory following task.
        In this case, the right hand side is a function from \f$\mathbf{R}\f$
        to \f$\mathbf{L}\f$.

        To set or get a time varying right hand side, use methods
        Implicit::rightHandSideFunction, Implicit::rightHandSideFromConfig.
        To get the value of the right hand
        side for a given time, use method Implicit::rightHandSideAt.
    */
    class HPP_CONSTRAINTS_DLLAPI Implicit {
      public:
	/// Operator equality
        bool operator== (const Implicit& other) const;
        /// Copy object and return shared pointer to copy
        virtual ImplicitPtr_t copy () const;
        /// Create a shared pointer to a new instance.
        /// \sa constructors
        static ImplicitPtr_t create
          (const DifferentiableFunctionPtr_t& func, ComparisonTypes_t comp,
	   std::vector<bool> mask = std::vector<bool>());

	/// Create a copy and return shared pointer
	static ImplicitPtr_t createCopy (const ImplicitPtr_t& other);

        virtual ~Implicit () {};

        /// \name Right hand side
        /// \{

        /// Computes the right hand side from a configuration
        ///
        /// in such a way that the configuration satisfies the numerical
        /// constraints
        /// \param config the input configuration.
        /// \retval rhs right hand side as a Lie group element.
        /// \precond rhs should be initialized with the right LiegroupSpace
        ///          (ie the output space of the function).
        ///
        /// \warning Only values of the right hand side corresponding to
        /// \link Equality "equality constraints" \endlink are set. As a
        /// result, the input configuration may not satisfy the other
        /// constraints.
        void rightHandSideFromConfig (ConfigurationIn_t config,
                                      LiegroupElementRef rhs);

        /// Check right hand side with regard to comparison types
        /// \param rhs  right hand side,
        /// \return true if right hand side is correct
        ///
        /// This method checks that elements of the right hand side log
        /// corresponding to comparison types different from Equality are equal
        /// to 0.
        bool checkRightHandSide(LiegroupElementConstRef rhs) const;

        /// Set the right hand side from a configuration
        ///
        /// in such a way that the configuration satisfies the numerical
        /// constraints
        /// \param config the input configuration.
        ///
        /// \warning Only values of the right hand side corresponding to
        /// \link Equality "equality constraints" \endlink are set. As a
        /// result, the input configuration may not satisfy the other
        /// constraints.
        ///
        /// \deprecated storing the right hand side in the constraint is
        /// deprecated. Users should handle themselves the right hand sides.
        /// Note that the solvers provided in this package already do.
        void rightHandSideFromConfig (ConfigurationIn_t config)
          HPP_CONSTRAINTS_DEPRECATED;

        /// Set the right hand side of the equation.
        /// \param rhs the right hand side.
        /// \deprecated storing the right hand side in the constraint is
        /// deprecated. Users should handle themselves the right hand sides.
        /// Note that the solvers provided in this package already do.
        void rightHandSide (vectorIn_t rhs) HPP_CONSTRAINTS_DEPRECATED;

        /// Return the right hand side of the equation.
        /// \deprecated storing the right hand side in the constraint is
        /// deprecated. Users should handle themselves the right hand sides.
        /// Note that the solvers provided in this package already do.
        vectorIn_t rightHandSide () const HPP_CONSTRAINTS_DEPRECATED;

        /// Get size of parameter defining the right hand side of the constraint
        ///
        /// See class documentation for details.
        size_type parameterSize () const;

        /// Get size of right hand side of the constraint
        ///
        /// This is the dimension (nq) of the output space of the function
        size_type rightHandSideSize () const;

        /// Return a modifiable reference to right hand side of equation.
        /// \deprecated storing the right hand side in the constraint is
        /// deprecated. Users should handle themselves the right hand sides.
        /// Note that the solvers provided in this package already do.
        vectorOut_t rightHandSide () HPP_CONSTRAINTS_DEPRECATED;

        /// Set time-varying right hand side
        /// \param rhsF Mapping from \f$\mathbf{R}\f$ to output space
        /// \f$\mathbf{L}\f$ of \f$h\f$ defining the variation along time
        /// of the right hand side.
        void rightHandSideFunction (const DifferentiableFunctionPtr_t& rhsF);

        /// Get time-varying right hand side
        /// \return the mapping from \f$\mathbf{R}\f$ to output space
        /// \f$\mathbf{L}\f$ of \f$h\f$ defining the variation along time
        /// of the right hand side.
        const DifferentiableFunctionPtr_t& rightHandSideFunction () const
        {
          return rhsFunction_;
        }

        /// Evaluate and set right hand side at given time
        /// \param s time
        vectorIn_t rightHandSideAt (const value_type& s);

        /// \}

        /// Return the ComparisonType
        const ComparisonTypes_t& comparisonType () const;

        /// Set the comparison type
        void comparisonType (const ComparisonTypes_t& comp);

        const segments_t& activeRows() const
        {
          return activeRows_;
        }

        /// Get indices of constraint coordinates that are equality
        const Eigen::RowBlockIndices& equalityIndices () const
        {
          return equalityIndices_;
        }

        /// Set inactive components of error to 0
        /// \retval error error vector computed by substracting a right hand
        ///         side to the output to the value of the function.
        void setInactiveRowsToZero(vectorOut_t error) const;

        /// Return a reference to function \f$h\f$.
        DifferentiableFunction& function () const
        {
          return *function_;
        }

        /// Return a reference to function \f$h\f$.
        const DifferentiableFunctionPtr_t& functionPtr () const
        {
          return function_;
        }

      protected:
        /// Constructor
        /// \param function the differentiable function
        /// \param comp vector of comparison \f$\mathbf{c}\f$.
        /// \param mask mask defining which components of the error are
        ///        taken into account to determine whether the constraint
        ///        is satisfied.
        /// \precond sizes of comp and of mask should be equal to size of
        ///          tangent space to
        ///          function output space \f$\mathbf{L}\f$.
        Implicit (const DifferentiableFunctionPtr_t& function,
                  ComparisonTypes_t comp, std::vector<bool> mask);

        /// Copy constructor
        Implicit (const Implicit& other);

        /// Test equality with other instance
        /// \param other object to copy
        /// \param swapAndTest whether we should also check other == this
        virtual bool isEqual (const Implicit& other, bool swapAndTest) const;

        // Store weak pointer to itself
        void init (const ImplicitWkPtr_t& weak)
        {
          weak_ = weak;
        }

        friend class ImplicitConstraintSet;
      private:
        void computeIndices();
        void computeActiveRows();
        ComparisonTypes_t comparison_;
        vector_t rhs_;
        size_type parameterSize_;
        DifferentiableFunctionPtr_t function_;
        DifferentiableFunctionPtr_t rhsFunction_;
        std::vector<bool> mask_;
        segments_t activeRows_;
        Eigen::RowBlockIndices inactiveRows_;
        std::vector<std::size_t> inequalityIndices_;
        Eigen::RowBlockIndices equalityIndices_;
        ImplicitWkPtr_t weak_;
        // To avoid dynamic memory allocation
        mutable LiegroupElement output_;
        mutable vector_t logOutput_;
      protected:
        Implicit () {}
      private:
        HPP_SERIALIZABLE();
    }; // class Implicit
    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_IMPLICIT_HH
