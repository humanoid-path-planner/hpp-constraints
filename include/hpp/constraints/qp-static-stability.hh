// Copyright (c) 2015, LAAS-CNRS
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
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

#ifndef HPP_CONSTRAINTS_QP_STATIC_STABILITY_HH
# define HPP_CONSTRAINTS_QP_STATIC_STABILITY_HH

# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/convex-shape-contact.hh>
# include <hpp/constraints/static-stability.hh>
# include <hpp/constraints/tools.hh>

# include "hpp/constraints/fwd.hh"

# include <qpOASES.hpp>

namespace hpp {
  namespace constraints {

    /// \addtogroup constraints
    /// \{

    class QPStaticStability : public DifferentiableFunction {
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

      private:
        qpOASES::real_t* Zeros;
        static const Eigen::Matrix <value_type, 6, 1> MinusGravity;
        const qpOASES::int_t nWSR;

        void impl_compute (vectorOut_t result, ConfigurationIn_t argument) const;

        void impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const;

        bool hasSolution (vectorOut_t dist) const;

        qpOASES::returnValue solveQP (vectorOut_t result) const;

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

        mutable qpOASES::SQProblem qp_;
        mutable MoE_t phi_;
        mutable qpOASES::real_t* A_;
        mutable InvertStorageOrderMap_t Amap_;
        mutable vector_t primal_, dual_;
    };
    /// \}
  } // namespace constraints
} // namespace hpp

#endif //  HPP_CONSTRAINTS_STATIC_STABILITY_HH
