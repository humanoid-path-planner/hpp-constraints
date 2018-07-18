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
    /// \addtogroup constraints
    /// \{

    /// This class represents a numerical constraint  with the following format
    /// \f{equation} f(q) comp rhs \f},
    /// where \f$comp\f$ is one of the following comparison operators:
    /// \li Equality: \f$f(\mathbf{q}) = rhs\f$, where \f$rhs\f$ is a
    /// parameterizable right hand side,
    /// \li EqualToZero: \f$f(\mathbf{q}) = 0\f$,
    /// \li Superior: \f$f(\mathbf{q}) > 0\f$
    /// \li Inferior: \f$f(\mathbf{q}) < 0\f$
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
        static ImplicitPtr_t create (const DifferentiableFunctionPtr_t& function);

        /// Create a shared pointer to a new instance.
        /// \sa constructors
        static ImplicitPtr_t create (const DifferentiableFunctionPtr_t& function,
            ComparisonTypes_t comp);

        /// Create a shared pointer to a new instance.
        /// \sa constructors
        static ImplicitPtr_t create (const DifferentiableFunctionPtr_t& function,
            ComparisonTypes_t comp, vectorIn_t rhs);

	/// Create a copy and return shared pointer
	static ImplicitPtr_t createCopy
	  (const ImplicitPtr_t& other);

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
        /// result, the input configuration may not satisfy the other constraints.
        /// The rationale is the following. Equality constraints define a
        /// foliation of the configuration space. Leaves of the foliation are
        /// defined by the value of the right hand side of the equality
        /// constraints. This method is mainly used in manipulation planning
        /// to retrieve the leaf a configuration lies on.
        void rightHandSideFromConfig (ConfigurationIn_t config);

        /// Set the right hand side of the equation.
        /// \param rhs the right hand side.
        void rightHandSide (vectorIn_t rhs);

        /// Return the right hand side of the equation.
        vectorIn_t rightHandSide () const;

        /// Return the size of the right hand side constraint.
        size_type rhsSize () const;

        /// Return the ComparisonType
        const ComparisonTypes_t& comparisonType () const;

	/// Set the comparison type
	void comparisonType (const ComparisonTypes_t& comp);

        bool constantRightHandSide () const;

        /// Return the right hand side of the equation.
        vectorOut_t nonConstRightHandSide ();

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
        vector_t& value ()
        {
          return value_;
        }

        /// Return a reference to the jacobian.
        /// This matrix can be used to store the derivative of the function,
        /// its size being initialized.
        matrix_t& jacobian ()
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

      private:
        ComparisonTypes_t comparison_;
        vector_t rhs_;
        size_type rhsRealSize_;
        DifferentiableFunctionPtr_t function_;
        vector_t value_;
        matrix_t jacobian_;
	ImplicitWkPtr_t weak_;
    }; // class Implicit
    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_IMPLICIT_HH
