// Copyright (c) 2014, LAAS-CNRS
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

#include "hpp/constraints/static-stability.hh"

#include <hpp/model/fcl-to-eigen.hh>

#include "hpp/constraints/orientation.hh"
#include <limits>
#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include "hpp/constraints/tools.hh"

namespace hpp {
  namespace constraints {
    namespace {
      static void cross (const fcl::Vec3f& v, eigen::matrix3_t& m)
      {
        m (0,1) = -v [2]; m (1,0) =  v [2];
        m (0,2) =  v [1]; m (2,0) = -v [1];
        m (1,2) = -v [0]; m (2,1) =  v [0];
        m (0,0) = m (1,1) = m (2,2) = 0;
      }
    }

    StaticStabilityGravity::StaticStabilityGravity ( const std::string& name,
        const DevicePtr_t& robot, const JointPtr_t& joint):
      DifferentiableFunction (robot->configSize (), robot->numberDof (), 5, name),
      robot_ (robot), joint_ (joint), jacobian_ (3, robot->numberDof ())
    {
      jacobian_.setZero ();
    }

    StaticStabilityGravityPtr_t StaticStabilityGravity::create (
        const std::string& name,
        const DevicePtr_t& robot,
        const JointPtr_t& joint)
    {
      return StaticStabilityGravityPtr_t (new StaticStabilityGravity (name, robot, joint));
    }

    StaticStabilityGravityPtr_t StaticStabilityGravity::create (const DevicePtr_t& robot,
        const JointPtr_t& joint)
    {
      return create ("StaticStabilityGravity", robot, joint);
    }

    void StaticStabilityGravity::addObjectTriangle (const fcl::TriangleP& t)
    {
      objectTriangles_.push_back (Triangle(t, joint_));
    }

    void StaticStabilityGravity::addFloorTriangle (const fcl::TriangleP& t)
    {
      floorTriangles_.push_back (Triangle(t));
    }

    void StaticStabilityGravity::impl_compute (vectorOut_t result, ConfigurationIn_t argument) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();

      selectTriangles ();

      fcl::Vec3f p = floor_->transformation ().transform (object_->center ());
      result [0] = p[0];
      if (floor_->isInside (floor_->intersection (object_->center (), floor_->normal ()))) {
        result [1] = 0;
        result [2] = 0;
      } else {
        result [1] = p[1];
        result [2] = p[2];
      }

      fcl::Matrix3f fTranspose = floor_->transformation ().getRotation ();
      fTranspose.transpose ();
      Rerror_ = object_->transformation().getRotation () * fTranspose;
      double theta;
      computeLog (r_, theta, Rerror_);
      result [3] = r_[1];
      result [4] = r_[2];
    }

    void StaticStabilityGravity::impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();

      selectTriangles ();

      const Transform3f& M = object_->joint_->currentTransformation ();
      const JointJacobian_t& Jjoint (object_->joint_->jacobian ());
      cross (M.getRotation () * object_->C_, Rcx_);
      eigen::matrix3_t Reigen;
      jacobian_.leftCols (Jjoint.cols ()) = floor_->transformation ().getRotation ()
        * ( - Rcx_ * Jjoint.bottomRows (3) + Jjoint.topRows (3));
      jacobian.row (0) = jacobian_.row (0);
      if (floor_->isInside (floor_->intersection (object_->center (), floor_->normal ()))) {
        jacobian.row (1).setZero ();
        jacobian.row (2).setZero ();
      } else {
        jacobian.row (1) = jacobian_.row (1);
        jacobian.row (2) = jacobian_.row (2);
      }

      fcl::Matrix3f fTranspose = floor_->transformation ().getRotation ();
      fTranspose.transpose ();
      Rerror_ = object_->transformation().getRotation () * fTranspose;
      double theta;
      computeLog (r_, theta, Rerror_);
      if (theta < 1e-3) {
	Jlog_.setIdentity ();
      } else {
        computeJlog (theta, r_, Jlog_);
      }
      jacobian_.leftCols (Jjoint.cols ()) = -Jlog_
        * object_->transformation().getRotation () * Jjoint.bottomRows (3);
      jacobian.row (3) = jacobian_.row (1);
      jacobian.row (4) = jacobian_.row (2);
    }

    void StaticStabilityGravity::selectTriangles () const
    {
      value_type dist, minDist = + std::numeric_limits <value_type>::infinity();
      for (Triangles::const_iterator o_it = objectTriangles_.begin ();
          o_it != objectTriangles_.end (); ++o_it) {
        o_it->updateToCurrentTransform ();
        const fcl::Vec3f& globalOC_ = o_it->center ();
        for (Triangles::const_iterator f_it = floorTriangles_.begin ();
            f_it != floorTriangles_.end (); ++f_it) {
          f_it->updateToCurrentTransform ();
          value_type dp = f_it->distance (f_it->intersection (globalOC_, f_it->normal ())),
                     dn = f_it->normal ().dot (globalOC_ - f_it->center ());
          if (dp < 0)
            dist = dn * dn;
          else
            dist = dp*dp + dn * dn;
          if (dist < minDist) {
            minDist = dist;
            object_ = o_it;
            floor_ = f_it;
          }
        }
      }
    }

    const value_type StaticStability::G = 9.81;
    const Eigen::Matrix <value_type, 6, 1> StaticStability::Gravity
      = (Eigen::Matrix <value_type, 6, 1>() << 0,0,-1, 0, 0, 0).finished();

    StaticStability::StaticStability ( const std::string& name,
        const DevicePtr_t& robot, const Contacts_t& contacts,
        const CenterOfMassComputationPtr_t& com):
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
          (1 + 6) * contacts.size() + 6, name),
      robot_ (robot), contacts_ (contacts), com_ (com),
      phi_ (Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero (6,contacts.size()),
          Eigen::Matrix<value_type, 6, Eigen::Dynamic>::Zero (6,contacts.size()*robot->numberDof()))
    {
      phi_.setSize (2,contacts.size());
      PointCom OG (com);
      for (std::size_t i = 0; i < contacts.size(); ++i) {
        PointInJoint OP1 (contacts[i].joint1,contacts[i].point1,robot->numberDof());
        PointInJoint OP2 (contacts[i].joint2,contacts[i].point2,robot->numberDof());
        VectorInJoint n1 (contacts[i].joint1,contacts[i].normal1,robot->numberDof()); 
        VectorInJoint n2 (contacts[i].joint2,contacts[i].normal2,robot->numberDof()); 

        phi_ (0,i) = CalculusBaseAbstract<>::create (n1);
        phi_ (1,i) = CalculusBaseAbstract<>::create ((OG - OP1) ^ n1);
        p1mp2s_.push_back (CalculusBaseAbstract<>::create (OP1 - OP2));
        n1mn2s_.push_back (CalculusBaseAbstract<>::create (n1 - n2));
      }
    }

    StaticStabilityPtr_t StaticStability::create ( const std::string& name,
        const DevicePtr_t& robot, const Contacts_t& contacts,
        const CenterOfMassComputationPtr_t& com)
    {
      return StaticStabilityPtr_t (new StaticStability (name, robot, contacts, com));
    }

    StaticStabilityPtr_t StaticStability::create (const DevicePtr_t& robot,
        const Contacts_t& contacts,
        const CenterOfMassComputationPtr_t& com)
    {
      return create ("StaticStability", robot, contacts, com);
    }

    void StaticStability::impl_compute (vectorOut_t result, ConfigurationIn_t argument) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();

      phi_.invalidate ();
      for (std::size_t i = 0; i < p1mp2s_.size(); ++i) {
        p1mp2s_[i]->invalidate();
        n1mn2s_[i]->invalidate();
      }
      phi_.computeValue (); //done in computePseudoInverse()
      phi_.computePseudoInverse ();
      const Eigen::Matrix <value_type, 6, 1> G = - 1 * Gravity;
      vector_t u = phi_.pinv() * G;
      // TODO: v should not be that but (V2 * V2*)^-1 u-
      vector_t v = 1 * (u.array () >= 0).select (0, -u);
      result.segment (0, contacts_.size()) = u + v - phi_.pinv() * (phi_.value() * v);
      result.segment <6> (contacts_.size()) = Gravity + phi_.value() * u;
      size_t shift = 6 + contacts_.size();
      for (std::size_t i = 0; i < p1mp2s_.size(); ++i) {
        p1mp2s_[i]->computeValue();
        n1mn2s_[i]->computeValue();
        result.segment <3> (shift    ) = p1mp2s_[i]->value();
        result.segment <3> (shift + 3) = n1mn2s_[i]->value();
        shift += 6;
      }
    }

    void StaticStability::impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();

      phi_.invalidate ();
      for (std::size_t i = 0; i < p1mp2s_.size(); ++i) {
        p1mp2s_[i]->invalidate();
        n1mn2s_[i]->invalidate();
      }
      phi_.computeValue (); //done in computePseudoInverse()
      phi_.computeJacobian ();
      phi_.computePseudoInverse ();

      const Eigen::Matrix <value_type, 6, 1> G = - 1 * Gravity;
      vector_t u = phi_.pinv() * G;
      matrix_t S = - matrix_t::Identity (u.size(), u.size());
      S.diagonal () = 1 * (u.array () >= 0).select
        (0, - vector_t::Ones (u.size()));
      vector_t v = S * u;

      Eigen::Matrix <value_type, 6, Eigen::Dynamic>
        JphiTimesV (6,robot_->numberDof());
      phi_.jacobianTimes (v, JphiTimesV);
      phi_.computePseudoInverseJacobian (G);
      jacobian.block (0, 0, contacts_.size(), robot_->numberDof()).noalias () =
        ( matrix_t::Identity (u.size(), u.size()) + S
          - phi_.pinv () * (phi_.value () * S)
          ) * phi_.pinvJacobian();
      jacobian.block (0, 0, contacts_.size(), robot_->numberDof()).noalias ()
        -= phi_.pinv () * JphiTimesV;
      phi_.computePseudoInverseJacobian (phi_.value () * v);
      jacobian.block (0, 0, contacts_.size(), robot_->numberDof()).noalias ()
        -= phi_.pinvJacobian ();

      phi_.jacobianTimes (u,
          jacobian.block (contacts_.size(), 0, 6, robot_->numberDof()));
      phi_.computePseudoInverseJacobian (Gravity);
      jacobian.block (contacts_.size(), 0, 6, robot_->numberDof())
        += - phi_.value() * phi_.pinvJacobian ();

      size_t shift = 6 + contacts_.size();
      for (std::size_t i = 0; i < p1mp2s_.size(); ++i) {
        p1mp2s_[i]->computeJacobian();
        n1mn2s_[i]->computeJacobian();
        jacobian.middleRows <3> (shift + i * 6    ) = p1mp2s_[i]->jacobian();
        jacobian.middleRows <3> (shift + i * 6 + 3) = n1mn2s_[i]->jacobian();
      }
    }
  } // namespace constraints
} // namespace hpp
