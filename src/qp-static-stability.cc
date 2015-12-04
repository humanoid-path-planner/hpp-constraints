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
    namespace {
      std::size_t forceDatasToNbContacts (
          const std::vector<ConvexShapeContact::ForceData>& fds) {
        std::size_t nb = 0;
        for (std::vector<ConvexShapeContact::ForceData>::const_iterator
            it = fds.begin (); it != fds.end (); ++it)
          nb += it->points.size ();
        return nb;
      }
    }

    const Eigen::Matrix <value_type, 6, 1> QPStaticStability::Gravity
      = (Eigen::Matrix <value_type, 6, 1>() << 0,0,-1, 0, 0, 0).finished();
    const Eigen::Matrix <value_type, 6, 1> QPStaticStability::MinusGravity
      = (Eigen::Matrix <value_type, 6, 1>() << 0,0,+1, 0, 0, 0).finished();

    QPStaticStability::QPStaticStability ( const std::string& name,
        const DevicePtr_t& robot, const Contacts_t& contacts,
        const CenterOfMassComputationPtr_t& com):
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
          1, name),
      Zeros (new qpOASES::real_t [contacts.size()]), nWSR (40),
      robot_ (robot), nbContacts_ (contacts.size()),
      com_ (com), H_ (nbContacts_,nbContacts_), G_ (nbContacts_),
      qp_ (nbContacts_, qpOASES::HST_SEMIDEF),
      phi_ (Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero (6,nbContacts_),
          Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero (6,nbContacts_*robot->numberDof())),
      primal_ (vector_t::Zero (nbContacts_)), dual_ (vector_t::Zero (nbContacts_))
    {
      VectorMap_t zeros (Zeros, nbContacts_); zeros.setZero ();

      qpOASES::Options options;
      qp_.setOptions( options );

      qp_.setPrintLevel (qpOASES::PL_NONE);
      phi_.setSize (2,nbContacts_);
      Traits<PointCom>::Ptr_t OG = PointCom::create(com);
      for (std::size_t i = 0; i < contacts.size(); ++i) {
        Traits<PointInJoint>::Ptr_t OP2 = PointInJoint::create
          (contacts[i].joint2,contacts[i].point2,robot->numberDof());
        Traits<VectorInJoint>::Ptr_t n2 = VectorInJoint::create
          (contacts[i].joint2,contacts[i].normal2,robot->numberDof()); 

        phi_ (0,i) = n2;
        phi_ (1,i) = (OG - OP2) ^ n2;
      }
    }

    QPStaticStability::QPStaticStability ( const std::string& name,
        const DevicePtr_t& robot, const std::vector<ForceData>& contacts,
        const CenterOfMassComputationPtr_t& com):
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
          1, name),
      Zeros (new qpOASES::real_t [forceDatasToNbContacts (contacts)]), nWSR (40),
      robot_ (robot), nbContacts_ (forceDatasToNbContacts (contacts)),
      com_ (com), H_ (nbContacts_, nbContacts_), G_ (nbContacts_),
      qp_ (nbContacts_, qpOASES::HST_SEMIDEF),
      phi_ (Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero (6,nbContacts_),
          Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero (6,nbContacts_*robot->numberDof())),
      primal_ (vector_t::Zero (nbContacts_)), dual_ (vector_t::Zero (nbContacts_))
    {
      VectorMap_t zeros (Zeros, nbContacts_); zeros.setZero ();

      qpOASES::Options options;
      qp_.setOptions( options );

      qp_.setPrintLevel (qpOASES::PL_NONE);
      phi_.setSize (2,nbContacts_);
      Traits<PointCom>::Ptr_t OG = PointCom::create(com);
      std::size_t col = 0;
      for (std::size_t i = 0; i < contacts.size (); ++i) {
        Traits<VectorInJoint>::Ptr_t n = VectorInJoint::create
          (contacts[i].joint,contacts[i].normal,robot->numberDof()); 
        for (std::size_t j = 0; j < contacts[i].points.size (); ++j) {
          Traits<PointInJoint>::Ptr_t OP = PointInJoint::create
            (contacts[i].joint,contacts[i].points[j],robot->numberDof());

          phi_ (0,col) = n;
          phi_ (1,col) = (OG - OP) ^ n;
          col++;
        }
      }
    }

    QPStaticStabilityPtr_t QPStaticStability::create ( const std::string& name,
        const DevicePtr_t& robot, const Contacts_t& contacts,
        const CenterOfMassComputationPtr_t& com)
    {
      return QPStaticStabilityPtr_t (new QPStaticStability (name, robot, contacts, com));
    }

    QPStaticStabilityPtr_t QPStaticStability::create ( const std::string& name,
        const DevicePtr_t& robot, const std::vector <ForceData>& contacts,
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

      qpOASES::returnValue ret = solveQP (result);
      if (ret != qpOASES::SUCCESSFUL_RETURN) {
        hppDout (error, "QP could not be solved. Error is " << ret);
      }
      if (!checkQPSol ()) {
        hppDout (error, "QP solution does not satisfies the constraints");
      }
    }

    void QPStaticStability::impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();

      phi_.invalidate ();
      phi_.computeSVD ();
      phi_.computeJacobian ();

      vector_t res(1);
      qpOASES::returnValue ret = solveQP (res);
      if (ret != qpOASES::SUCCESSFUL_RETURN) {
        hppDout (error, "QP could not be solved. Error is " << ret);
      }
      if (!checkQPSol ()) {
        hppDout (error, "QP solution does not satisfies the constraints");
      }
      if (!checkStrictComplementarity ()) {
        hppDout (error, "Strict complementary slackness does not hold. "
            "Jacobian WILL be wrong.");
      }

      matrix_t JT_phi_F (nbContacts_, robot_->numberDof());
      matrix_t J_F (6, robot_->numberDof());
      phi_.jacobianTransposeTimes (phi_.value () * primal_, JT_phi_F);
      phi_.jacobianTimes (primal_, J_F);

      jacobian = 0.5 * primal_.transpose() * JT_phi_F
        + (0.5 * phi_.value() * primal_ + Gravity).transpose() * J_F;
    }

    inline qpOASES::returnValue QPStaticStability::solveQP
      (vectorOut_t result) const
    {
      // TODO: Use the SVD to solve a smaller quadratic problem
      // Try to find a positive solution
      using qpOASES::SUCCESSFUL_RETURN;

      H_ = phi_.value ().transpose () * phi_.value();
      G_ = phi_.value().transpose () * Gravity;

      qpOASES::int_t nwsr = nWSR;
      qp_.reset ();
      qp_.setHessianType (qpOASES::HST_SEMIDEF);
      qpOASES::returnValue ret;
      ret = qp_.init (H_.data(), G_.data(), Zeros, 0, nwsr, 0);
      qp_.getPrimalSolution (primal_.data ());
      qp_.getDualSolution (dual_.data ());
      result[0] = 2*qp_.getObjVal () + MinusGravity.squaredNorm ();
      return ret;
    }

    bool QPStaticStability::checkQPSol () const
    {
      return (primal_.array () >= -1e-8).all();
    }

    bool QPStaticStability::checkStrictComplementarity () const
    {
      qpOASES::real_t eps = qp_.getOptions().boundTolerance;
      return (
             (primal_.array() > eps &&   dual_.cwiseAbs().array() <= eps )
          || (  dual_.array() > eps && primal_.cwiseAbs().array() <= eps )
          ).all();
    }
  } // namespace constraints
} // namespace hpp
