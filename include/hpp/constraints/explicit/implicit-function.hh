// Copyright (c) 2015 - 2018 LAAS-CNRS
// Authors: Florent Lamiraux, Joseph Mirabel
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
    ///  f, g are differentiable functions with values in  a Lie group.
    ///
    ///  This class is mainly used to create hpp::constraints::Explicit
    ///  instances.
    class ImplicitFunction : public DifferentiableFunction
    {
    public:
      /// create instance and return shared pointer
      /// \deprecated used create method that takes a LiegroupSpace instead
      ///             of a robot as input.
      typedef boost::shared_ptr <ImplicitFunction> Ptr_t;
      static Ptr_t create
      (const DevicePtr_t& robot, const DifferentiableFunctionPtr_t& function,
       const segments_t& inputConf, const segments_t& outputConf,
       const segments_t& inputVelocity, const segments_t& outputVelocity)
        HPP_CONSTRAINTS_DEPRECATED;

      /// create instance and return shared pointer
      /// \deprecated used create method that takes a LiegroupSpace instead
      ///             of a robot as input.
      static Ptr_t create
      (const DevicePtr_t& robot, const DifferentiableFunctionPtr_t& function,
       const DifferentiableFunctionPtr_t& g,
       const segments_t& inputConf, const segments_t& outputConf,
       const segments_t& inputVelocity, const segments_t& outputVelocity)
        HPP_CONSTRAINTS_DEPRECATED;

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
      /// \deprecated used constructor that takes a LiegroupSpace instead
      ///             of a robot as input.
      ImplicitFunction (const DevicePtr_t& robot,
			const DifferentiableFunctionPtr_t& function,
			const segments_t& inputConf,
                        const segments_t& outputConf,
			const segments_t& inputVelocity,
			const segments_t& outputVelocity)
        HPP_CONSTRAINTS_DEPRECATED;

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
    }; // class ImplicitFunction

    } // namespace explicit_
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_EXPLICIT_IMPLICIT_FUNCTION_HH
