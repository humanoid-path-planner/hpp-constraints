//
// Copyright (c) 2014 CNRS
// Authors: Florent Lamiraux
//
// This file is part of hpp-constraints
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
// hpp-constraints  If not, see
// <http://www.gnu.org/licenses/>.

#ifndef HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_HH
# define HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_HH

# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/config.hh>
# include <hpp/pinocchio/liegroup-element.hh>

namespace hpp {
  namespace constraints {

    /// \addtogroup constraints
    /// \{

    /// Differentiable function
    class HPP_CONSTRAINTS_DLLAPI DifferentiableFunction
    {
    public:
      virtual ~DifferentiableFunction () {}
      /// Evaluate the function at a given parameter.
      ///
      /// \note parameters should be of the correct size.
      LiegroupElement operator () (vectorIn_t argument) const
      {
	assert (argument.size () == inputSize ());
        LiegroupElement result (outputSpace_);
	impl_compute (result, argument);
        return result;
      }
      /// Evaluate the function at a given parameter.
      ///
      /// \note parameters should be of the correct size.
      void value (LiegroupElementRef result,
                  vectorIn_t argument) const
      {
	assert (result.size () == outputSize ());
	assert (argument.size () == inputSize ());
	impl_compute (result, argument);
      }
      /// Computes the jacobian.
      ///
      /// \retval jacobian jacobian will be stored in this argument
      /// \param argument point at which the jacobian will be computed
      void jacobian (matrixOut_t jacobian, vectorIn_t argument) const
      {
	assert (argument.size () == inputSize ());
	assert (jacobian.rows () == outputDerivativeSize ());
	assert (jacobian.cols () == inputDerivativeSize ());
	impl_jacobian (jacobian, argument);
      }

      /// Returns a vector of booleans that indicates whether the corresponding
      /// configuration parameter influences this constraints.
      const ArrayXb& activeParameters () const
      {
        return activeParameters_;
      }

      /// Returns a vector of booleans that indicates whether the corresponding
      /// velocity parameter influences this constraints.
      const ArrayXb& activeDerivativeParameters () const
      {
        return activeDerivativeParameters_;
      }

      /// Get dimension of input vector
      size_type inputSize () const
      {
	return inputSize_;
      }
      /// Get dimension of input derivative vector
      ///
      /// The dimension of configuration vectors might differ from the dimension
      /// of velocity vectors since some joints are represented by non minimal
      /// size vectors: e.g. quaternion for SO(3)
      size_type inputDerivativeSize () const
      {
	return inputDerivativeSize_;
      }
      /// Get output element
      ///
      /// Getting an output element enables users to know the type of Liegroup
      /// the function output values lie in.
      LiegroupSpacePtr_t outputSpace () const
      {
        return outputSpace_;
      }
      /// Get dimension of output vector
      size_type  outputSize () const
      {
	return outputSpace_->nq ();
      }
      /// Get dimension of output derivative vector
      size_type  outputDerivativeSize () const
      {
	return outputSpace_->nv ();
      }
      /// \brief Get function name.
      ///
      /// \return Function name.
      const std::string& name () const
      {
	return name_;
      }

      /// Display object in a stream
      virtual std::ostream& print (std::ostream& o) const;

      std::string context () const {
        return context_;
      }

      void context (const std::string& c) {
        context_ = c;
      }

      /// Approximate the jacobian using forward finite difference.
      /// \retval jacobian jacobian will be stored in this argument
      /// \param arg point at which the jacobian will be computed
      /// \param robot use to add configuration and velocities. If set to NULL,
      ///              the configuration space is considered a vector space.
      /// \param eps refers to \f$\epsilon\f$ in
      ///            http://en.wikipedia.org/wiki/Numerical_differentiation
      /// Evaluate the function (x.size() + 1) times but less precise the
      /// finiteDifferenceCentral
      void finiteDifferenceForward (matrixOut_t jacobian, vectorIn_t arg,
          DevicePtr_t robot = DevicePtr_t (),
          value_type eps = std::sqrt(Eigen::NumTraits<value_type>::epsilon())) const;

      /// Approximate the jacobian using forward finite difference.
      /// \retval jacobian jacobian will be stored in this argument
      /// \param arg point at which the jacobian will be computed
      /// \param robot use to add configuration and velocities. If set to NULL,
      ///              the configuration space is considered a vector space.
      /// \param eps refers to \f$\epsilon\f$ in
      ///            http://en.wikipedia.org/wiki/Numerical_differentiation
      /// Evaluate the function 2*x.size() times but more precise the
      /// finiteDifferenceForward
      void finiteDifferenceCentral (matrixOut_t jacobian, vectorIn_t arg,
          DevicePtr_t robot = DevicePtr_t (),
          value_type eps = std::sqrt(Eigen::NumTraits<value_type>::epsilon())) const;

    protected:
      /// \brief Concrete class constructor should call this constructor.
      ///
      /// \param sizeInput dimension of the function input
      /// \param sizeInputDerivative dimension of the function input derivative,
      /// \param sizeOutput dimension of the output,
      /// \param name function's name
      DifferentiableFunction (size_type sizeInput,
			      size_type sizeInputDerivative,
			      size_type sizeOutput,
			      std::string name = std::string ());

      /// \brief Concrete class constructor should call this constructor.
      ///
      /// \param sizeInput dimension of the function input
      /// \param sizeInputDerivative dimension of the function input derivative,
      /// \param outputSpace output space of the function.
      /// \param name function name
      DifferentiableFunction (size_type sizeInput,
			      size_type sizeInputDerivative,
			      const LiegroupSpacePtr_t& outputSpace,
			      std::string name = std::string ());

      /// User implementation of function evaluation
      virtual void impl_compute (LiegroupElementRef result,
				 vectorIn_t argument) const = 0;

      virtual void impl_jacobian (matrixOut_t jacobian,
				  vectorIn_t arg) const = 0;

      /// Dimension of input vector.
      size_type inputSize_;
      /// Dimension of input derivative
      size_type inputDerivativeSize_;
      /// Dimension of output vector
      LiegroupSpacePtr_t outputSpace_;

      /// Initialized to true by this class. Child class are responsible for
      /// updating it.
      /// \sa activeParameters
      ArrayXb activeParameters_;
      /// Initialized to true by this class. Child class are responsible for
      /// updating it.
      /// \sa activeDerivativeParameters
      ArrayXb activeDerivativeParameters_;

    private:
      std::string name_;
      /// Context of creation of function
      std::string context_;

      friend class DifferentiableFunctionSet;
    }; // class DifferentiableFunction
    inline std::ostream&
    operator<< (std::ostream& os, const DifferentiableFunction& f)
    {
      return f.print (os);
    }
    /// \}
  } // namespace constraints
} // namespace hpp


#endif // HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_HH
