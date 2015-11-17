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

      Eigen::Matrix
        <qpOASES::real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        Hessian (const std::size_t nbC)
      {
        typedef Eigen::Matrix<qpOASES::real_t, Eigen::Dynamic, Eigen::Dynamic,
                Eigen::RowMajor> RowMajorMatrix_t;
        RowMajorMatrix_t ret = RowMajorMatrix_t::Zero (2*nbC,2*nbC);
        ret.block (nbC, nbC, nbC, nbC).setIdentity ();
        ret.block (0, 0, nbC, nbC).diagonal ().setConstant (1e-6);
        return ret;
      }
    }

    QPStaticStabilityFeasabilityPtr_t QPStaticStabilityFeasability::create
      (QPStaticStabilityPtr_t stab)
    {
      return QPStaticStabilityFeasabilityPtr_t (
          new QPStaticStabilityFeasability (stab));
    }

    QPStaticStabilityFeasability::QPStaticStabilityFeasability
      (QPStaticStabilityPtr_t stab) :
        DifferentiableFunction (stab->inputSize (), stab->inputDerivativeSize(),
            6, stab->name () + "_stability"),
        stab_ (stab)
        {}

    void QPStaticStabilityFeasability::impl_compute
      (vectorOut_t result, ConfigurationIn_t argument) const
    {
      stab_->robot_->currentConfiguration (argument);
      stab_->robot_->computeForwardKinematics ();

      stab_->phi_.invalidate ();
      stab_->phi_.computeSVD ();
      stab_->hasSolution (result);
    }

    void QPStaticStabilityFeasability::impl_jacobian
      (matrixOut_t jacobian, ConfigurationIn_t argument) const
    {
      stab_->robot_->currentConfiguration (argument);
      stab_->robot_->computeForwardKinematics ();

      stab_->phi_.invalidate ();
      stab_->phi_.computeSVD ();
      stab_->phi_.computeJacobian ();

      vector_t sol = stab_->phi_.svd().solve (QPStaticStability::MinusGravity);
      stab_->phi_.jacobianTimes (sol, jacobian);
      stab_->phi_.computePseudoInverseJacobian (QPStaticStability::MinusGravity);
      jacobian.noalias () += stab_->phi_.value() * stab_->phi_.pinvJacobian ();
    }

    const Eigen::Matrix <value_type, 6, 1> QPStaticStability::Gravity
      = (Eigen::Matrix <value_type, 6, 1>() << 0,0,-1, 0, 0, 0).finished();
    const Eigen::Matrix <value_type, 6, 1> QPStaticStability::MinusGravity
      = (Eigen::Matrix <value_type, 6, 1>() << 0,0,+1, 0, 0, 0).finished();

    QPStaticStability::QPStaticStability ( const std::string& name,
        const DevicePtr_t& robot, const Contacts_t& contacts,
        const CenterOfMassComputationPtr_t& com):
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
          contacts.size(), name),
      Zeros (new qpOASES::real_t [2*contacts.size()]), nWSR (40),
      robot_ (robot), nbContacts_ (contacts.size()),
      com_ (com), H_ (Hessian (nbContacts_)),
      qp_ (2*nbContacts_, 6, qpOASES::HST_IDENTITY),
      phi_ (Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero (6,nbContacts_),
          Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero (6,nbContacts_*robot->numberDof())),
      A_ (new qpOASES::real_t [6*2*nbContacts_]), Amap_ (A_, 6, 2*nbContacts_),
      primal_ (vector_t::Zero (2*nbContacts_)), dual_ (vector_t::Zero (6 + 2*nbContacts_))
    {
      VectorMap_t zeros (Zeros, 2*nbContacts_); zeros.setZero ();

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
          forceDatasToNbContacts (contacts), name),
      Zeros (new qpOASES::real_t [2*forceDatasToNbContacts (contacts)]), nWSR (20),
      robot_ (robot), nbContacts_ (forceDatasToNbContacts (contacts)),
      com_ (com), H_ (Hessian (nbContacts_)),
      qp_ (2*nbContacts_, 6, qpOASES::HST_IDENTITY),
      phi_ (Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero (6,nbContacts_),
          Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero (6,nbContacts_*robot->numberDof())),
      A_ (new qpOASES::real_t [6*2*nbContacts_]), Amap_ (A_, 6, 2*nbContacts_),
      primal_ (vector_t::Zero (2*nbContacts_)), dual_ (vector_t::Zero (6 + 2*nbContacts_))
    {
      VectorMap_t zeros (Zeros, 2*nbContacts_); zeros.setZero ();

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

      if (solveQP (result) != qpOASES::SUCCESSFUL_RETURN) {
        if (!checkQPSol ()) {
          hppDout (error, "QPStaticStability should not be used when "
              "QPStaticStabilityFeasability is not satisfied");
          result.setConstant (-1);
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

      vector_t sol (nbContacts_);
      if (solveQP (sol) != qpOASES::SUCCESSFUL_RETURN) {
        hppDout (error, "QPStaticStability should not be used when "
            "QPStaticStabilityFeasability is not satisfied");
        jacobian.setZero ();
        return;
      }
      matrix_t jtf (6,robot_->numberDof());
      matrix_t tmp (6,robot_->numberDof());
      phi_.jacobianTimes (primal_.segment (0, nbContacts_), jtf);
      phi_.jacobianTimes (primal_.segment (nbContacts_, nbContacts_), tmp);
      jtf.noalias() -= tmp;
      vector_t diag = vector_t::Ones(2*nbContacts_);

      // phi_.jacobianTransposeTimes (dual_.segment <6> (2*nbContacts_),
          // jacobian);
      // jacobian *= -1;
      qpOASES::Bounds b;
      qp_.getBounds (b);
      const qpOASES::Indexlist* il = b.getFixed ();
      for (qpOASES::int_t i = 0; i < il->getLength (); ++i) {
        diag (il->getNumber (i)) = 0;
        // if (il->getNumber (i) >= nbContacts_)
          // jacobian.row (il->getNumber (i) - nbContacts_).setZero ();
      }

      matrix_t pi (Amap_.cols (), Amap_.rows ());
      MoE_t::SVD_t svd (Amap_ * diag.asDiagonal(), Eigen::ComputeThinU | Eigen::ComputeThinV);
      // svd.setThreshold (1e-4);
      pseudoInverse <MoE_t::SVD_t> (svd, pi);
      matrix_t pk (pi.rows(),pi.rows());
      projectorOnKernel <MoE_t::SVD_t> (svd, pk);
      jacobian = - pk.bottomRows (nbContacts_) * H_ * pi * jtf;
    }

    inline bool QPStaticStability::hasSolution (vectorOut_t dist) const
    {
      dist.noalias () = getU2 <MoE_t::SVD_t> (phi_.svd()) *
        ( getU2 <MoE_t::SVD_t> (phi_.svd()).adjoint() * Gravity );
      return dist.squaredNorm () < 1e-8;
    }

    inline qpOASES::returnValue QPStaticStability::solveQP
      (vectorOut_t result) const
    {
      // TODO: Use the SVD to solve a smaller quadratic problem
      // Try to find a positive solution
      // using qpOASES::QProblem;
      using qpOASES::HST_IDENTITY;
      using qpOASES::SUCCESSFUL_RETURN;

      Amap_.leftCols  (nbContacts_) =  phi_.value (); // Need a copy because of the row-major order.
      Amap_.rightCols (nbContacts_) = -phi_.value (); // Need a copy because of the row-major order.

      qpOASES::int_t nwsr = nWSR;
      const qpOASES::real_t eps = 1e-4;
      const qpOASES::real_t lbA[6] = {-eps,-eps,1-eps,-eps,-eps,-eps}; // - StaticStability::Gravity
      const qpOASES::real_t ubA[6] = { eps, eps,1+eps, eps, eps, eps}; // - StaticStability::Gravity
      qp_.reset ();
      qpOASES::returnValue ret;
      if (qp_.isInitialised()) {
        ret =
          qp_.hotstart (H_.data(), Zeros, A_, Zeros, 0, lbA, ubA, nwsr, 0);
      } else {
        ret =
          qp_.init (H_.data(), Zeros, A_, Zeros, 0, lbA, ubA, nwsr, 0, primal_.data());
      }
      qp_.getPrimalSolution (primal_.data ());
      qp_.getDualSolution (dual_.data ());
      result = primal_.segment (nbContacts_, nbContacts_);
      return ret;
    }

    bool QPStaticStability::checkQPSol () const
    {
      vector_t error = Amap_ * primal_;
      error.noalias () += Gravity;
      return error.isZero (6*1e-4);
    }
  } // namespace constraints
} // namespace hpp
