//
// Copyright (c) 2014 CNRS
// Authors: Florent Lamiraux
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

#ifndef HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_HH
# define HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_HH

#include <hpp/util/debug.hh>
# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/config.hh>
# include <hpp/pinocchio/liegroup-element.hh>

# include <hpp/util/serialization-fwd.hh>

namespace hpp {
  namespace constraints {

    /// \addtogroup constraints
    /// \{

    /// Differentiable function from a \link hpp::pinocchio::LiegroupSpace Lie
    /// group\endlink, for instance the configuration space of a robot
    /// (hpp::pinocchio::Device) to a another \link
    /// hpp::pinocchio::LiegroupSpace Lie group\endlink.
    ///
    /// Note that the input Lie group is only represented by the sizes
    /// of the elements and of the velocities: methods \link
    /// DifferentiableFunction::inputSize inputSize\endlink and \link
    /// DifferentiableFunction::inputDerivativeSize
    /// inputDerivativeSize\endlink
    ///
    /// The output space can be accessed by method \link
    /// DifferentiableFunction::outputSpace outputSpace\endlink.
    ///
    /// The value of the function for a given input can be accessed by method
    /// \link DifferentiableFunction::value value \endlink.
    /// The Jacobian of the function for a given input can be accessed by
    /// method \link DifferentiableFunction::jacobian jacobian \endlink.
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
	assert (result.space()->nq() == outputSize ());
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
      /// Get output space
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

      bool operator== (DifferentiableFunction const & other) const;
      bool operator!= (DifferentiableFunction const & b) const;

      /// Return pair of joints the relative pose between which
      /// modifies the function value if any
      ///
      /// This method is useful to know whether an instance of Implicit constrains
      /// the relative pose between two joints.
      /// \return the pair of joints involved, arranged in order of increasing
      /// joint index, or a pair of empty shared pointers.
      /// \note
      ///   \li if absolute pose (relative pose with respect to "universe"),
      ///       "universe" is returned as empty shared pointer
      ///   \li child class reimplementing this may require a valid "robot"
      ///       argument, which the constraints are applied on.
      virtual std::pair<JointConstPtr_t, JointConstPtr_t> dependsOnRelPoseBetween
          (DeviceConstPtr_t /*robot*/) const
      {
        return std::pair<JointConstPtr_t, JointConstPtr_t>(nullptr, nullptr);
      };

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

      virtual bool isEqual(const DifferentiableFunction& other) const
      {
        if (name_ != other.name_)
          return false;
        if (inputSize_ != other.inputSize_)
          return false;
        if (inputDerivativeSize_ != other.inputDerivativeSize_)
          return false;
        if (*outputSpace_ != *(other.outputSpace_))
          return false;

        return true;
      }

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

    protected:
      DifferentiableFunction() {}
    private:
      HPP_SERIALIZABLE();
    }; // class DifferentiableFunction
    inline std::ostream&
    operator<< (std::ostream& os, const DifferentiableFunction& f)
    {
      return f.print (os);
    }
    /// \}
  } // namespace constraints
} // namespace hpp

/*
// This will override boost::shared_ptr 's equality operator
// between 2 DifferentiableFunctionPtr_t
namespace boost {
  using namespace hpp::constraints;
  typedef DifferentiableFunction T;
  typedef shared_ptr<T> TPtr;

  template<> inline bool operator==<T, T> (TPtr const & a, TPtr const & b) BOOST_SP_NOEXCEPT
  {
    return *a == *b;
  }

  template<> inline bool operator!=<T, T> (TPtr const & a, TPtr const & b) BOOST_SP_NOEXCEPT
  {
    return !(a == b);
  }
}
*/

#endif // HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_HH
