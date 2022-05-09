// Copyright (c) 2015 - 2018 LAAS-CNRS
// Authors: Florent Lamiraux, Joseph Mirabel
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

#ifndef HPP_CONSTRAINTS_EXPLICIT_IMPLICIT_FUNCTION_HH
# define HPP_CONSTRAINTS_EXPLICIT_IMPLICIT_FUNCTION_HH

# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/matrix-view.hh>

namespace hpp {
  namespace constraints {
    namespace explicit_ {
    /// Function of the form q -> g (q_out) - f (q_in)
    ///
    /// where
    ///  \li q_out is a vector composed of configuration variables of
    ///      q,
    ///  \li q_in is the vector composed of other configuration variables of
    ///      q,
    ///  \li f, g are differentiable functions with values in a Lie group.
    ///
    ///  This class is mainly used to create hpp::constraints::Explicit
    ///  instances.
    class ImplicitFunction : public DifferentiableFunction
    {
    public:
      typedef shared_ptr <ImplicitFunction> Ptr_t;
      /// create instance and return shared pointer
      ///
      /// \param configSpace input space of this function - usually a robot
      ///                    configuration space,
      /// \param function function f,
      /// \param inputConf set of indices defining q_in,
      /// \param inputVelocity set of indices defining q_in derivative,
      /// \param outputConf set of indices defining q_out
      /// \param outputVelocity set of indices defining q_out derivative
      ///

      static Ptr_t create
      (const LiegroupSpacePtr_t& configSpace,
       const DifferentiableFunctionPtr_t& function,
       const segments_t& inputConf, const segments_t& outputConf,
       const segments_t& inputVelocity, const segments_t& outputVelocity);

      /// Get function f that maps input variables to output variables
      const DifferentiableFunctionPtr_t& inputToOutput () const;

    protected:
      /// Constructor
      /// \param configSpace input space of this function - usually a robot
      ///                    configuration space,
      /// \param function function f,
      /// \param inputConf set of indices defining q_in,
      /// \param inputVelocity set of indices defining q_in derivative,
      /// \param outputConf set of indices defining q_out
      /// \param outputVelocity set of indices defining q_out derivative
      ImplicitFunction (const LiegroupSpacePtr_t& configSpace,
                const DifferentiableFunctionPtr_t& function,
                const segments_t& inputConf,
                const segments_t& outputConf,
                const segments_t& inputVelocity,
                const segments_t& outputVelocity);
      /// Compute g (q_out) - f (q_in)
      void impl_compute (LiegroupElementRef result, vectorIn_t argument) const;

      /// Compute Jacobian of g (q_out) - f (q_in) with respect to q.
      void impl_jacobian (matrixOut_t jacobian, vectorIn_t arg) const;

      bool isEqual(const DifferentiableFunction& other) const {
        const ImplicitFunction& castother = dynamic_cast<const ImplicitFunction&>(other);
        if (!DifferentiableFunction::isEqual(other))
          return false;

        if (robot_ != castother.robot_)
          return false;
        if (inputToOutput_ != castother.inputToOutput_)
          return false;
        if (inputConfIntervals_.rows() != castother.inputConfIntervals_.rows())
          return false;
        if (outputConfIntervals_.rows() != castother.outputConfIntervals_.rows())
          return false;
        if (inputDerivIntervals_.rows() != castother.inputDerivIntervals_.rows())
          return false;
        if (outputDerivIntervals_.rows() != castother.outputDerivIntervals_.rows())
          return false;

        return true;
      }

      /// Return pair of joints the relative pose between which
      /// modifies the function value if any
      ///
      /// Currently checks if the implicit function specifies a joint
      /// where
      ///  \li q_out is a vector corresponding to only 1 joint
      ///  \li q_in is an empty vector (since f is constant and specifies
      ///      the whole or part of the pose of the joint)
      ///
      /// \param robot the robot the constraints are applied on,
      /// \return pair of pointers to the lock joint and its parent joint,
      /// arranged in order of increasing joint index, or a pair of empty
      /// shared pointers if the implicit function does not specify a locked
      /// joint.
      std::pair<JointConstPtr_t, JointConstPtr_t> dependsOnRelPoseBetween
          (DeviceConstPtr_t robot) const;

    private:
      void computeJacobianBlocks ();

      DevicePtr_t robot_;
      DifferentiableFunctionPtr_t inputToOutput_;
      Eigen::RowBlockIndices inputConfIntervals_;
      Eigen::RowBlockIndices outputConfIntervals_;
      Eigen::RowBlockIndices inputDerivIntervals_;
      Eigen::RowBlockIndices outputDerivIntervals_;
      std::vector <Eigen::MatrixBlocks <false, false> > outJacobian_;
      std::vector <Eigen::MatrixBlocks <false, false> > inJacobian_;
      mutable vector_t qIn_;
      mutable LiegroupElement f_qIn_, qOut_;
      mutable LiegroupElement result_;
      // Jacobian of explicit function
      mutable matrix_t Jf_;

      ImplicitFunction() {}
      HPP_SERIALIZABLE();
    }; // class ImplicitFunction

    } // namespace explicit_
  } // namespace constraints
} // namespace hpp

BOOST_CLASS_EXPORT_KEY(hpp::constraints::explicit_::ImplicitFunction)

#endif // HPP_CONSTRAINTS_EXPLICIT_IMPLICIT_FUNCTION_HH
