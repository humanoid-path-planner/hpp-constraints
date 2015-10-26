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
#include <hpp/model/eigen.hh>
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
          contacts.size() + 6, name),
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

      qp_.setPrintLevel (qpOASES::PL_NONE);
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
      phi_.computeSVD ();

      // Is there a solution to the unconstrained problem ?
      bool hasSol = hasSolution (result.segment <6> (0));

      if (hasSol) {
        if (solveQP (result.segment (6, contacts_.size()))
            == qpOASES::SUCCESSFUL_RETURN)
          return;
        else {
          // There are no positive solution.
          result.segment (6, contacts_.size ()).noalias () =
            phi_.svd ().solve (MinusGravity);
        }
      } else {
        // No solution.
        // TODO: Does this introduce a discontinuity ?
        result.segment (6, contacts_.size ()).noalias () =
          phi_.svd ().solve (MinusGravity);
      }
    }

    void QPStaticStability::impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();

      phi_.invalidate ();
      phi_.computeSVD ();
      phi_.computeJacobian ();

      // Is there a solution to the unconstrained problem ?
      vector_t u (6);
      bool hasSol = hasSolution (u);

      vector_t sol = phi_.svd().solve (MinusGravity);
      phi_.jacobianTimes (sol, jacobian.block (0, 0, 6, robot_->numberDof()));
      phi_.computePseudoInverseJacobian (MinusGravity);
      jacobian.block (0, 0, 6, robot_->numberDof()).noalias ()
        += phi_.value() * phi_.pinvJacobian ();

      if (hasSol) {
        if (solveQP (sol) == qpOASES::SUCCESSFUL_RETURN) {
          phi_.jacobianTransposeTimes (dual_.segment <6> (contacts_.size()),
              jacobian.block (6, 0, contacts_.size(), robot_->numberDof()));
          qpOASES::Bounds b;
          qp_.getBounds (b);
          const qpOASES::Indexlist* il = b.getFixed ();
          for (qpOASES::int_t i = 0; i < il->getLength (); ++i)
            jacobian.row (6 + il->getNumber (i)).setZero ();
        }
        else {
          // There are no positive solution.
          // phi_.computePseudoInverseJacobian (MinusGravity);
          jacobian.block (6, 0, contacts_.size(), robot_->numberDof())
            = phi_.pinvJacobian ();
        }
      } else {
        // phi_.computePseudoInverseJacobian (MinusGravity);
        jacobian.block (6, 0, contacts_.size(), robot_->numberDof()) =
          phi_.pinvJacobian ();
      }
    }

    inline bool QPStaticStability::hasSolution (vectorOut_t dist) const
    {
      using namespace hpp::model;

      dist.noalias () = getU2 <MoE_t::SVD_t> (phi_.svd()) *
        ( getU2 <MoE_t::SVD_t> (phi_.svd()).adjoint() * Gravity );
      return dist.squaredNorm () < 1e-6;
    }

    inline qpOASES::returnValue QPStaticStability::solveQP
      (vectorOut_t result) const
    {
      // TODO: Use the SVD to solve a smaller quadratic problem
      // Try to find a positive solution
      // using qpOASES::QProblem;
      using qpOASES::HST_IDENTITY;
      using qpOASES::SUCCESSFUL_RETURN;

      Amap_ = phi_.value (); // Need a copy because of the row-major order.

      qpOASES::int_t nwsr = nWSR;
      const qpOASES::real_t eps = 1e-4;
      const qpOASES::real_t lbA[6] = {-eps,-eps,1-eps,-eps,-eps,-eps}; // - StaticStability::Gravity
      const qpOASES::real_t ubA[6] = { eps, eps,1+eps,-eps, eps, eps}; // - StaticStability::Gravity
      qp_.reset ();
      qp_.setHessianType (HST_IDENTITY);
      qpOASES::returnValue ret;
      if (qp_.isInitialised()) {
        ret =
          qp_.hotstart (NULL, Zeros, A_, Zeros, 0, lbA, ubA, nwsr, 0);
      } else {
        ret =
          qp_.init (NULL, Zeros, A_, Zeros, 0, lbA, ubA, nwsr, 0, primal_.data());
      }
      qp_.getPrimalSolution (primal_.data ());
      qp_.getDualSolution (dual_.data ());
      result = primal_;
      return ret;
    }
  } // namespace constraints
} // namespace hpp
