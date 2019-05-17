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

namespace hpp {
  namespace constraints {
    /**
        \addtogroup constraints
        \{
    */
    /**
        This class represents a parameterizable numerical constraint that
        compares the output of a function \f$h\f$ to a right hand side vector.

        \par Definition
        \li The function \f$h\f$ takes input in a configuration space
        \f$\mathcal{C}\f$ and output in a Lie group \f$\mathbf{L}\f$,
        \li the dimensions of \f$\mathbf{L}\f$ and of its tangent space are
        respectively \f$(n_q,n_v)\f$.
        \li The comparison is represented by a vector \f$c\f$ of dimension
        \f$n_v\f$
        with values in enum hpp::contraints::ComparisonType = {
        \f$\mathbf{Equality}\f$, \f$\mathbf{EqualToZero}\f$,
        \f$\mathbf{Inferior}\f$, \f$\mathbf{Superior}\f$ }.
        \li The right hand side is a vector of dimension \f$n_q\f$.

        \par Error
        A configuration \f$\mathbf{q}\f$ is said to satisfy the constraint for
        a given right hand side if and only if the error \f$e\f$ as computed
        below is
        smaller in norm than a given threshold.

        Let
        \f[
        \Delta = f (\mathbf{q}) - rhs \in \mathbf{R}^{n_v},
        \f]
        for each component \f$i\in\{0,\cdots,n_v-1\}\f$,
        \li if \f$c_i\f$ is \f$\mathbf{Inferior}\f$,
        \f$e_i = \max (0,\Delta_i)\f$,
        \li if \f$c_i\f$ is \f$\mathbf{Superior}\f$,
        \f$e_i = \min (0,\Delta_i)\f$,
        \li if \f$c_i\f$ is \f$\mathbf{Equality}\f$,
        \f$e_i = \Delta_i\f$,
        \li if \f$c_i\f$ is \f$\mathbf{EqualToZero}\f$, \f$e_i = \Delta_i\f$.

        \par Parameterizable right hand side
        Lines with \textbf{Equality} comparator in the above definition of the
        error need a parameter, while lines with comparators \textbf{Inferior},
        \textbf{Superior}, or \textbf{EqualToZero} do not.
        As a consequence, the right hand side of the constraint is defined by
        a vector \f$\lambda\f$ of parameters of size the number of
        \textbf{Equality} occurences in vector \f$c\f$. The right hand side is
        then defined as in the following example:
        \f[
        rhs = \exp\left(\begin{array}{c}\lambda_1 \\ 0 \\ 0 \\ \lambda_2 \\
        \vdots \end{array}\right)
         \ \ \ \ c = \left(\begin{array}{c}\mathbf{Equality} \\
        \mathbf{EqualToZero} \\ \mathbf{Inferior} \\ \mathbf{Equality} \\
        \vdots \end{array}\right)
        \f]
        To retrieve the size of vector \f$\lambda\f$, call method
        Implicit::parameterSize ().
    */
    class HPP_CONSTRAINTS_DLLAPI Implicit {
      public:
	/// Operator equality
	bool operator== (const Implicit& other) const
	{
	  return isEqual (other, true);
	}
        /// Copy object and return shared pointer to copy
        virtual ImplicitPtr_t copy () const;
        /// Create a shared pointer to a new instance.
        /// \sa constructors
        static ImplicitPtr_t create
          (const DifferentiableFunctionPtr_t& function);

        /// Create a shared pointer to a new instance.
        /// \sa constructors
        static ImplicitPtr_t create
          (const DifferentiableFunctionPtr_t& function, ComparisonTypes_t comp);

        /// Create a shared pointer to a new instance.
        /// \sa constructors
        static ImplicitPtr_t create
          (const DifferentiableFunctionPtr_t& function,
           ComparisonTypes_t comp, vectorIn_t rhs);

	/// Create a copy and return shared pointer
	static ImplicitPtr_t createCopy (const ImplicitPtr_t& other);

        virtual ~Implicit () {};

        /// \name Right hand side
        /// \{

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
        /// \deprecated In future versions, right hand side will not be a member
        ///             of the class anymore.
        void rightHandSideFromConfig (ConfigurationIn_t config) HPP_CONSTRAINTS_DEPRECATED;

        /// Set the right hand side of the equation.
        /// \param rhs the right hand side.
        ///
        /// \deprecated In future versions, right hand side will not be a member
        ///             of the class anymore.
        void rightHandSide (vectorIn_t rhs) HPP_CONSTRAINTS_DEPRECATED;

        /// Return the right hand side of the equation.
        ///
        /// \deprecated In future versions, right hand side will not be a member
        ///             of the class anymore.
        vectorIn_t rightHandSide () const HPP_CONSTRAINTS_DEPRECATED;

        /// Return the dimension of the left hand side function tangent output
        /// space
        ///
        /// \deprecated Use function ()->outputSpace ()->nv () instead.
        size_type rhsSize () const HPP_CONSTRAINTS_DEPRECATED;

        /// Get size of parameter defining the right hand side of the constraint
        ///
        /// See class documentation for details.
        size_type parameterSize () const;

        /// Get size of right hand side of the constraint
        ///
        /// This is the dimension (nq) of the output space of the function
        size_type rightHandSideSize () const;

        /// Return the ComparisonType
        const ComparisonTypes_t& comparisonType () const;

	/// Set the comparison type
	void comparisonType (const ComparisonTypes_t& comp);

        bool constantRightHandSide () const;

        /// Return a modifiable reference to right hand side of equation.
        /// \deprecated In future versions, right hand side will not be a member
        ///             of the class anymore.
        vectorOut_t nonConstRightHandSide () HPP_CONSTRAINTS_DEPRECATED;

        void rightHandSideFunction (const DifferentiableFunctionPtr_t& rhsF);

        const DifferentiableFunctionPtr_t& rightHandSideFunction () const
        {
          return rhsFunction_;
        }

        vectorIn_t rightHandSideAt (const value_type& s);

        /// \}

        /// Return a reference to function \f$f\f$.
        DifferentiableFunction& function () const
        {
          return *function_;
        }

        /// Return a reference to the inner function.
        const DifferentiableFunctionPtr_t& functionPtr () const
        {
          return function_;
        }

        /// Return a reference to the value.
        /// This vector can be used to store the output of the function,
        /// its size being initialized.
        /// \deprecated In future versions, the value will not be a member of
        ///             the class anymore.
        vector_t& value () HPP_CONSTRAINTS_DEPRECATED
        {
          return value_;
        }

        /// Return a reference to the jacobian.
        /// This matrix can be used to store the derivative of the function,
        /// its size being initialized.
        /// \deprecated In future versions, the Jacobian will not be a member of
        ///             the class anymore.
        matrix_t& jacobian () HPP_CONSTRAINTS_DEPRECATED
        {
          return jacobian_;
        }

      protected:
        /// Constructor
        /// \param function the differentiable function
        Implicit (const DifferentiableFunctionPtr_t& function,
            ComparisonTypes_t comp);

        /// Constructor
        /// \param function the differentiable function
        /// \param rhs the right hand side of this equation
        Implicit (const DifferentiableFunctionPtr_t& function,
            ComparisonTypes_t comp, vectorIn_t rhs);

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
        ComparisonTypes_t comparison_;
        vector_t rhs_;
        size_type rhsRealSize_;
        DifferentiableFunctionPtr_t function_;
        DifferentiableFunctionPtr_t rhsFunction_;
        vector_t value_;
        matrix_t jacobian_;
	ImplicitWkPtr_t weak_;
    }; // class Implicit
    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_IMPLICIT_HH
