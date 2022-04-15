// Copyright (c) 2014, LAAS-CNRS
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

#ifndef HPP_CONSTRAINTS_STATIC_STABILITY_HH
# define HPP_CONSTRAINTS_STATIC_STABILITY_HH

# include <vector>

# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/config.hh>
# include <hpp/constraints/deprecated.hh>

# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/symbolic-calculus.hh>

namespace hpp {
  namespace constraints {

    /// \addtogroup constraints
    /// \{
    class HPP_CONSTRAINTS_DLLAPI StaticStability :
      public DifferentiableFunction {
      public:
        static const value_type G;
        static const Eigen::Matrix <value_type, 6, 1> Gravity;

        struct Contact_t {
          JointPtr_t joint;
          vector3_t point;
          vector3_t normal;
          bool operator==(Contact_t const& other) const {
            if (joint != other.joint) return false;
            if (point != other.point) return false;
            if (normal != other.normal) return false;
            return true;
          }
        };
        typedef std::vector <Contact_t> Contacts_t;

        /// Constructor
        /// \param robot the robot the constraints is applied to,
        /// \param com COM of the object in the joint frame.
        StaticStability (const std::string& name, const DevicePtr_t& robot,
            const Contacts_t& contacts,
            const CenterOfMassComputationPtr_t& com);

        static StaticStabilityPtr_t create (
            const std::string& name,
            const DevicePtr_t& robot,
            const Contacts_t& contacts,
            const CenterOfMassComputationPtr_t& com);

        static StaticStabilityPtr_t create (
            const DevicePtr_t& robot,
            const Contacts_t& contacts,
            const CenterOfMassComputationPtr_t& com);

        MatrixOfExpressions<>& phi () {
          return phi_;
        }

      protected:
        bool isEqual(const DifferentiableFunction& other) const {
          const StaticStability& castother = dynamic_cast<const StaticStability&>(other);
          if (!DifferentiableFunction::isEqual(other))
            return false;
          
          if (robot_ != castother.robot_)
            return false;
          if (contacts_ != castother.contacts_)
            return false;
          if (com_ != castother.com_)
            return false;
          
          return true;
        }
      private:
        void impl_compute (LiegroupElementRef result,
                           ConfigurationIn_t argument) const;

        void impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const;

        static void findBoundIndex (vectorIn_t u, vectorIn_t v, 
            value_type& lambdaMin, size_type* iMin,
            value_type& lambdaMax, size_type* iMax);

        /// Return false if uMinus.isZero(), i which case v also zero (not computed).
        bool computeUminusAndV (vectorIn_t u, vectorOut_t uMinus,
            vectorOut_t v) const;

        void computeVDot (const ConfigurationIn_t arg, vectorIn_t uMinus, vectorIn_t S,
            matrixIn_t uDot, matrixOut_t uMinusDot, matrixOut_t vDot) const;

        void computeLambdaDot (vectorIn_t u, vectorIn_t v, const std::size_t i0,
            matrixIn_t uDot, matrixIn_t vDot, vectorOut_t lambdaDot) const;

        DevicePtr_t robot_;
        Contacts_t contacts_;
        CenterOfMassComputationPtr_t com_;

        typedef MatrixOfExpressions<eigen::vector3_t, JacobianMatrix> MoE_t;

        mutable MoE_t phi_;
        mutable vector_t u_, uMinus_, v_;
        mutable matrix_t uDot_, uMinusDot_, vDot_;
        mutable vector_t lambdaDot_; 
    };
    /// \}
  } // namespace constraints
} // namespace hpp

#endif //  HPP_CONSTRAINTS_STATIC_STABILITY_HH
