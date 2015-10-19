// Copyright (c) 2015, Joseph Mirabel
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

#include "hpp/constraints/qp-static-stability.hh"

#include <hpp/model/fcl-to-eigen.hh>

#include <limits>
#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include "hpp/constraints/tools.hh"

namespace hpp {
  namespace constraints {
    const Eigen::Matrix <value_type, 6, 1> QPStaticStability::Gravity
      = (Eigen::Matrix <value_type, 6, 1>() << 0,0,-1, 0, 0, 0).finished();
    const Eigen::Matrix <value_type, 6, 1> QPStaticStability::MinusGravity
      = (Eigen::Matrix <value_type, 6, 1>() << 0,0,+1, 0, 0, 0).finished();

    QPStaticStability::QPStaticStability ( const std::string& name,
        const DevicePtr_t& robot, const Contacts_t& contacts,
        const CenterOfMassComputationPtr_t& com):
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
          contacts.size(), name),
      Zeros (new qpOASES::real_t [contacts.size()]), nWSR (20),
      robot_ (robot), contacts_ (contacts), com_ (com),
      qp_ (contacts.size(), 6, qpOASES::HST_IDENTITY),
      phi_ (Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero (6,contacts.size()),
          Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero (6,contacts.size()*robot->numberDof())),
      A_ (new qpOASES::real_t [6*contacts.size()]), Amap_ (A_, 6, contacts.size()),
      primal_ (vector_t::Zero (contacts.size())), dual_ (vector_t::Zero (6 + contacts.size()))
    {
      VectorMap_t zeros (Zeros, contacts.size()); zeros.setZero ();

      qpOASES::Options options;
      qp_.setOptions( options );

      qp_.setPrintLevel (qpOASES::PL_LOW);
      phi_.setSize (2,contacts.size());
      PointCom OG (com);
      for (std::size_t i = 0; i < contacts.size(); ++i) {
        PointInJoint OP1 (contacts[i].joint1,contacts[i].point1,robot->numberDof());
        PointInJoint OP2 (contacts[i].joint2,contacts[i].point2,robot->numberDof());
        VectorInJoint n1 (contacts[i].joint1,contacts[i].normal1,robot->numberDof()); 
        VectorInJoint n2 (contacts[i].joint2,contacts[i].normal2,robot->numberDof()); 

        phi_ (0,i) = CalculusBaseAbstract<>::create (n2);
        phi_ (1,i) = CalculusBaseAbstract<>::create ((OG - OP2) ^ n2);
      }
    }

    QPStaticStabilityPtr_t QPStaticStability::create ( const std::string& name,
        const DevicePtr_t& robot, const Contacts_t& contacts,
        const CenterOfMassComputationPtr_t& com)
    {
      return QPStaticStabilityPtr_t (new QPStaticStability (name, robot, contacts, com));
    }

    QPStaticStabilityPtr_t QPStaticStability::create (const DevicePtr_t& robot,
        const Contacts_t& contacts,
        const CenterOfMassComputationPtr_t& com)
    {
      return create ("QPStaticStability", robot, contacts, com);
    }

    void QPStaticStability::impl_compute (vectorOut_t result, ConfigurationIn_t argument) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();

      phi_.invalidate ();
      phi_.computeValue ();

      Amap_ = phi_.value (); // Need a copy because of the row-major order.

      // const value_type g[8] = {0,0,0,0,0,0,0,0}; // - StaticStability::Gravity
      // const value_type lb[8] = {0,0,0,0,0,0,0,0}; // - StaticStability::Gravity
      qpOASES::int_t nwsr = nWSR;
      const qpOASES::real_t eps = 1e-4;
      const qpOASES::real_t lbA[6] = {-eps,-eps,1-eps,-eps,-eps,-eps}; // - StaticStability::Gravity
      const qpOASES::real_t ubA[6] = { eps, eps,1+eps,-eps, eps, eps}; // - StaticStability::Gravity

      // Assume problem is feasible
      qp_.reset ();
      qp_.setHessianType (qpOASES::HST_IDENTITY);
      if (qpOASES::SUCCESSFUL_RETURN == qp_.init (NULL, Zeros, A_, Zeros, 0,
          lbA, ubA, nwsr, 0, primal_.data())) {
        // , dual_.data());

        qp_.getPrimalSolution (primal_.data());
        qp_.getDualSolution (dual_.data());
        result.segment (0, contacts_.size()) = primal_;
      } else {
        qp_.reset ();
        qp_.setHessianType (qpOASES::HST_IDENTITY);
        /// Same problem with no bound on F
        if (qpOASES::SUCCESSFUL_RETURN == qp_.init (NULL, Zeros, A_, 0, 0,
              lbA, ubA, nwsr, 0, primal_.data())) {
          qp_.getPrimalSolution (primal_.data());
          qp_.getDualSolution (dual_.data());
        }
      }
    }

    void QPStaticStability::impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();

      phi_.invalidate ();

      phi_.computeSVD ();
      phi_.computeJacobian ();

      // Assume problem is feasible
      // TODO: Do not need to solve the QP again.
      qp_.reset ();
      qp_.setHessianType (qpOASES::HST_IDENTITY);
      Amap_ = phi_.value (); // Need a copy because of the row-major order.
      qpOASES::int_t nwsr = nWSR;
      const qpOASES::real_t eps = 1e-4;
      const qpOASES::real_t lbA[6] = {-eps,-eps,1-eps,-eps,-eps,-eps}; // - StaticStability::Gravity
      const qpOASES::real_t ubA[6] = { eps, eps,1+eps,-eps, eps, eps}; // - StaticStability::Gravity
      qp_.init (NULL, Zeros, A_, Zeros, 0,
          lbA, ubA, nwsr, 0, primal_.data());
          // , dual_.data());
      qp_.getPrimalSolution (primal_.data());
      qp_.getDualSolution (dual_.data());

      phi_.jacobianTransposeTimes (dual_.segment <6> (contacts_.size()),
          jacobian.block (0, 0, contacts_.size(), robot_->numberDof()));
    }
  } // namespace constraints
} // namespace hpp
