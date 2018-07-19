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

#ifndef HPP_CONSTRAINTS_EXPLICIT_FUNCTION_HH
# define HPP_CONSTRAINTS_EXPLICIT_FUNCTION_HH

# include <hpp/constraints/differentiable-function.hh>

namespace hpp {
  namespace constraints {
    namespace explicit_ {
    /// Function of the form f (q) = q2 - g (q1)
    ///
    /// where
    ///  \li q2 is a vector composed of a subset of configuration variables of
    ///      q,
    ///  \li q1 is the vector composed of the other configuration variables of
    ///      q,
    ///  g is a differentiable function with values in  a Lie group.
    ///
    ///  This class is mainly used to create hpp::constraints::Explicit
    ///  instances.
    template <bool GisIdentity>
    class Function : public DifferentiableFunction
    {
    public:
      typedef boost::shared_ptr <Function> Ptr_t;
      static Ptr_t create
      (const DevicePtr_t& robot, const DifferentiableFunctionPtr_t& function,
       const segments_t& inputConf, const segments_t& inputVelocity,
       const segments_t& outputConf, const segments_t& outputVelocity);

      static Ptr_t create
      (const DevicePtr_t& robot, const DifferentiableFunctionPtr_t& function,
       const DifferentiableFunctionPtr_t& g,
       const segments_t& inputConf, const segments_t& inputVelocity,
       const segments_t& outputConf, const segments_t& outputVelocity);

      /// Get function that maps input variables to output variables
      const DifferentiableFunctionPtr_t& inputToOutput () const;

    protected:
      Function (const DevicePtr_t& robot,
			const DifferentiableFunctionPtr_t& function,
                        const DifferentiableFunctionPtr_t& g,
			const segments_t& inputConf,
			const segments_t& inputVelocity,
                        const segments_t& outputConf,
			const segments_t& outputVelocity);

      /// Compute q_{output} - f (q_{input})
      void impl_compute (LiegroupElement& result, vectorIn_t argument) const;

      /// Compute Jacobian of q_{output} - f (q_{input})
      void impl_jacobian (matrixOut_t jacobian, vectorIn_t arg) const;

    private:
      void computeJacobianBlocks ();

      struct GenericGData {
        DifferentiableFunctionPtr_t g_;
        mutable LiegroupElement g_qOut_;
        mutable matrix_t Jg_;

        GenericGData (const DifferentiableFunctionPtr_t& g)
          : g_ (g), g_qOut_ (g->outputSpace()),
          Jg_ (g->outputSpace()->nv(), g->inputDerivativeSize())
        {}
        void computeValue (const LiegroupElement& qOut) const
        {
          g_->value (g_qOut_, qOut.vector());
        }
        const LiegroupElement& value (const LiegroupElement&) const
        {
          return g_qOut_;
        }
        matrix_t& jacobian (const LiegroupElement& qOut) const
        {
          g_->jacobian (Jg_, qOut.vector());
          return Jg_;
        }
      };

      struct IdentityData {
        mutable matrix_t Jg_;
        IdentityData (const DifferentiableFunctionPtr_t&) {}
        void computeValue (const LiegroupElement&) const {}
        const LiegroupElement& value (const LiegroupElement& qOut) const { return qOut; }
        const matrix_t& jacobian (const LiegroupElement&) const { return Jg_; }
      };

      typedef typename boost::conditional<GisIdentity, IdentityData, GenericGData>::type GData;

      DevicePtr_t robot_;
      DifferentiableFunctionPtr_t inputToOutput_;
      Eigen::RowBlockIndices inputConfIntervals_;
      Eigen::RowBlockIndices inputDerivIntervals_;
      Eigen::RowBlockIndices outputConfIntervals_;
      Eigen::RowBlockIndices outputDerivIntervals_;
      std::vector <Eigen::MatrixBlocks <false, false> > outJacobian_;
      std::vector <Eigen::MatrixBlocks <false, false> > inJacobian_;
      GData gData_;
      mutable vector_t qIn_;
      mutable LiegroupElement f_qIn_, qOut_;
      mutable LiegroupElement result_;
      // Jacobian of explicit function
      mutable matrix_t Jf_;
    }; // class Function

    typedef Function<true > BasicFunction;
    typedef Function<false> GenericFunction;

  } // namespace explicit_
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_EXPLICIT_FUNCTION_HH
