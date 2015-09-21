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
						     const DevicePtr_t& robot):
      DifferentiableFunction (robot->configSize (), robot->numberDof (), 5,
			      name),
      robot_ (robot), jacobian_ (3, robot->numberDof ())
    {
      jacobian_.setZero ();
    }

    StaticStabilityGravityPtr_t StaticStabilityGravity::create (
        const std::string& name,
        const DevicePtr_t& robot)
    {
      return StaticStabilityGravityPtr_t (new StaticStabilityGravity
					  (name, robot));
    }

    StaticStabilityGravityPtr_t StaticStabilityGravity::create
    (const DevicePtr_t& robot)
    {
      return create ("StaticStabilityGravity", robot);
    }

    void StaticStabilityGravity::addObjectTriangle (const fcl::TriangleP& t,
						    const JointPtr_t& joint)
    {
      objectTriangles_.push_back (std::make_pair (Triangle(t), joint));
    }

    void StaticStabilityGravity::addFloorTriangle (const fcl::TriangleP& t,
						   const JointPtr_t& joint)
    {
      floorTriangles_.push_back (std::make_pair (Triangle(t), joint));
    }

    void StaticStabilityGravity::impl_compute (vectorOut_t result, ConfigurationIn_t argument) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();

      selectTriangles ();

      const Transform3f M (inverse (floor_->second->currentTransformation ()) *
			   object_->second->currentTransformation ());
      fcl::Vec3f p = floor_->first->inversePosition ().transform
	(M.transform (object_->center ()));
      result [0] = p[0];
      if (floor_->isInside (floor_->intersection (M.transform (object_->center ()), floor_->normal ()))) {
        result [1] = 0;
        result [2] = 0;
      } else {
        result [1] = p[1];
        result [2] = p[2];
      }

      Rerror_ = object_->second->inversePosition ().getRotation () *
        transpose ((floor_->inversePosition ().getRotation () *
		    M.getRotation ()));
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

      const Transform3f M (inverse (floor_->second->currentTransformation ()) *
			   object_->second->currentTransformation ());
      const JointJacobian_t& Jo (object_->second->jacobian ());
      cross (M.getRotation () * object_->center (), Rcx_);
      eigen::matrix3_t RTeigen;
      do {
        const fcl::Matrix3f& R = floor_->inversePosition ().getRotation ();
        RTeigen << R (0,0), R (0, 1), R (0, 2),
                   R (1,0), R (1, 1), R (1, 2),
                   R (2,0), R (2, 1), R (2, 2);
      } while (0);
      jacobian_.leftCols (Jo.cols ()) = RTeigen *
        (- Rcx_ * Jo.bottomRows (3) + Jo.topRows (3));
      jacobian.row (0) = jacobian_.row (0);
      if (floor_->isInside (floor_->intersection (M.transform (object_->center ()), floor_->normal ()))) {
        jacobian.row (1).setZero ();
        jacobian.row (2).setZero ();
      } else {
        jacobian.row (1) = jacobian_.row (1);
        jacobian.row (2) = jacobian_.row (2);
      }
      Rerror_ = object_->inversePosition ().getRotation () * (floor_->inversePosition ().getRotation () * M.getRotation ()).transpose ();
      double theta;
      computeLog (r_, theta, Rerror_);
      if (theta < 1e-3) {
	Jlog_.setIdentity ();
      } else {
        computeJlog (theta, r_, Jlog_);
      }
      do {
        fcl::Matrix3f RT = M.getRotation ();
        RT.transpose ();
        const fcl::Matrix3f R = object_->inversePosition ().getRotation () * RT;
        RTeigen << R (0,0), R (0, 1), R (0, 2),
                   R (1,0), R (1, 1), R (1, 2),
                   R (2,0), R (2, 1), R (2, 2);
      } while (0);
      jacobian_.leftCols (Jo.cols ()) = -Jlog_ * RTeigen * Jo.bottomRows (3);
      jacobian.row (3) = jacobian_.row (1);
      jacobian.row (4) = jacobian_.row (2);
    }

    void StaticStabilityGravity::selectTriangles () const
    {
      fcl::Vec3f globalOC_, globalFC;

      value_type dist, minDist = + std::numeric_limits <value_type>::infinity();
      for (Triangles::const_iterator o_it = objectTriangles_.begin ();
          o_it != objectTriangles_.end (); ++o_it) {
	const Transform3f& Mo = o_it->second->currentTransformation ();
        globalOC_ = Mo.transform (o_it->first->center ());
        for (Triangles::const_iterator f_it = floorTriangles_.begin ();
            f_it != floorTriangles_.end (); ++f_it) {
	  const Transform3f& Mf = f_it->second->currentTransformation ();
	  const Matrix3_t& Rf = Mf->getRotation ();
	  vector3_t nf (Rf * f_it->first->normal ());
	  globalFC = Mf * f_it->first->center ();
          value_type dp = f_it->first->distance
	    (f_it->first->intersection (globalOC_, nf)),
	    dn = nf.dot (globalOC_ - globalFC);
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
      result.segment (0, contacts_.size()) = - /*com_->mass() * G * */ phi_.pinv() * Gravity;
      result.segment <6> (contacts_.size()) = /* G * */ Gravity - phi_.value() * (phi_.pinv() * /* G * */ Gravity);
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
      phi_.computePseudoInverseJacobian (/* G * */ Gravity);
      jacobian.block (0, 0, contacts_.size(), robot_->numberDof()) =
        - /* com_->mass() * */ phi_.pinvJacobian();
      phi_.jacobianTimes (- phi_.pinv() * /* G * */ Gravity,
          jacobian.block (contacts_.size(), 0, 6, robot_->numberDof()));
      jacobian.block (contacts_.size(), 0, 6, robot_->numberDof())
        += - phi_.value() * phi_.pinvJacobian ();
      //jacobian.block (contacts_.size(), 0, contacts_.size(), robot_->numberDof()) =
        //- phi_.jacobianTimes (phi_.pinv() * Gravity)
        //- phi_.value() * phi_.pinvJacobian ();
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
