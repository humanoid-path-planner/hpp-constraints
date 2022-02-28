// Copyright (c) 2015, LAAS-CNRS
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

#ifndef HPP_CONSTRAINTS_QP_STATIC_STABILITY_HH
# define HPP_CONSTRAINTS_QP_STATIC_STABILITY_HH

# include <hpp/constraints/fwd.hh>

# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/convex-shape-contact.hh>
# include <hpp/constraints/static-stability.hh>

# include <qpOASES.hpp>

namespace hpp {
  namespace constraints {

    /// \addtogroup constraints
    /// \{
    class HPP_CONSTRAINTS_DLLAPI QPStaticStability : public DifferentiableFunction {
      public:
        static const Eigen::Matrix <value_type, 6, 1> Gravity;

        typedef StaticStability::Contact_t Contact_t;
        typedef StaticStability::Contacts_t Contacts_t;
        typedef ConvexShapeContact::ForceData ForceData;

        /// Constructor
        /// \param robot the robot the constraints is applied to,
        /// \param com COM of the robot
        QPStaticStability (const std::string& name, const DevicePtr_t& robot,
            const Contacts_t& contacts,
            const CenterOfMassComputationPtr_t& com);

        /// Constructor
        /// \param robot the robot the constraints is applied to,
        /// \param com COM of the robot
        QPStaticStability (const std::string& name, const DevicePtr_t& robot,
            const std::vector <ForceData>& contacts,
            const CenterOfMassComputationPtr_t& com);

        static QPStaticStabilityPtr_t create (
            const std::string& name,
            const DevicePtr_t& robot,
            const Contacts_t& contacts,
            const CenterOfMassComputationPtr_t& com);

        static QPStaticStabilityPtr_t create (
            const std::string& name,
            const DevicePtr_t& robot,
            const std::vector <ForceData>& contacts,
            const CenterOfMassComputationPtr_t& com);

        static QPStaticStabilityPtr_t create (
            const DevicePtr_t& robot,
            const Contacts_t& contacts,
            const CenterOfMassComputationPtr_t& com);

        MatrixOfExpressions<>& phi () {
          return phi_;
        }
      protected:
        bool isEqual(const DifferentiableFunction& other) const {
          const QPStaticStability& castother = dynamic_cast<const QPStaticStability&>(other);
          if (!DifferentiableFunction::isEqual(other))
            return false;
          
          if (name() != castother.name())
            return false;
          if (robot_ != castother.robot_)
            return false;
          if (nbContacts_ != castother.nbContacts_)
            return false;
          if (com_ != castother.com_)
            return false;
          if (Zeros != castother.Zeros)
            return false;
          if (nWSR != castother.nWSR)
            return false;
          
          return true;
        }
      private:
        static const Eigen::Matrix <value_type, 6, 1> MinusGravity;

        qpOASES::real_t* Zeros;
        const qpOASES::int_t nWSR;

        void impl_compute (LiegroupElementRef result, ConfigurationIn_t argument)
          const;

        void impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const;

        qpOASES::returnValue solveQP (vectorOut_t result) const;

        bool checkQPSol () const;
        bool checkStrictComplementarity () const;

        DevicePtr_t robot_;
        std::size_t nbContacts_;
        CenterOfMassComputationPtr_t com_;

        typedef MatrixOfExpressions<eigen::vector3_t, JacobianMatrix> MoE_t;
        typedef Eigen::Matrix<qpOASES::real_t, Eigen::Dynamic, Eigen::Dynamic,
                Eigen::RowMajor> RowMajorMatrix_t;
        typedef Eigen::Map <RowMajorMatrix_t> InvertStorageOrderMap_t;
        typedef Eigen::Map <Eigen::Matrix <qpOASES::real_t, Eigen::Dynamic, 1>
          > VectorMap_t;
        typedef Eigen::Map <const vector_t> ConstVectorMap_t;

        mutable RowMajorMatrix_t H_;
        mutable vector_t G_;
        mutable qpOASES::QProblemB qp_;
        mutable MoE_t phi_;
        mutable vector_t primal_, dual_;
    };
    /// \}
  } // namespace constraints
} // namespace hpp

#endif //  HPP_CONSTRAINTS_STATIC_STABILITY_HH
