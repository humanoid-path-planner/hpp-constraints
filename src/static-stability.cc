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
#include "hpp/constraints/relative-transformation.hh"
#include <limits>
#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include "hpp/constraints/tools.hh"

namespace hpp {
  namespace constraints {

    StaticStabilityGravity::StaticStabilityGravity
    (const std::string& name, const DevicePtr_t& robot) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (), 5,
			      name),
      robot_ (robot), relativeTransformation_ (RelativeTransformation::create
					       (name, robot,
						robot->rootJoint (),
						robot->rootJoint (),
						Transform3f (), Transform3f (),
						boost::assign::list_of (true)
						(true)(true)(true)(true)(true))
					       )
    {
      result_.resize (6);
      jacobian_.resize (6, robot->numberDof ());
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
      (*relativeTransformation_) (result_, argument);
      result [0] = result_ [0];
      result [1] = result_ [1];
      result [2] = result_ [2];
      result [3] = result_ [4];
      result [4] = result_ [5];
      if (isInside_) {
	result [1] = 0;
        result [2] = 0;
      }
      hppDout (info, "result = " << result.transpose ());
    }

    void StaticStabilityGravity::computeInternalJacobian
    (ConfigurationIn_t argument) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      selectTriangles ();
      relativeTransformation_->jacobian (jacobian_, argument);
    }

    void StaticStabilityGravity::impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const
    {
      computeInternalJacobian (argument);
      jacobian.row (0) = jacobian_.row (0);
      jacobian.row (1) = jacobian_.row (1);
      jacobian.row (2) = jacobian_.row (2);
      jacobian.row (3) = jacobian_.row (4);
      jacobian.row (4) = jacobian_.row (5);
      if (isInside_) {
	jacobian.row (0) = jacobian_.row (0);
	jacobian.row (1).setZero ();
	jacobian.row (2).setZero ();
	jacobian.row (3) = jacobian_.row (4);
	jacobian.row (4) = jacobian_.row (5);
      } else {
	jacobian.row (0) = jacobian_.row (0);
	jacobian.row (1) = jacobian_.row (1);
	jacobian.row (2) = jacobian_.row (2);
	jacobian.row (3) = jacobian_.row (4);
	jacobian.row (4) = jacobian_.row (5);
      }
    }

    void StaticStabilityGravity::selectTriangles () const
    {
      fcl::Vec3f globalOC_;
      Triangles::const_iterator object;
      Triangles::const_iterator floor;

      value_type dist, minDist = + std::numeric_limits <value_type>::infinity();
      for (Triangles::const_iterator o_it = objectTriangles_.begin ();
          o_it != objectTriangles_.end (); ++o_it) {
        for (Triangles::const_iterator f_it = floorTriangles_.begin ();
            f_it != floorTriangles_.end (); ++f_it) {
	  Transform3f Mf;
	  if (f_it->second) {
	    Mf = f_it->second->currentTransformation ();
	  }
	  Transform3f M (inverse (Mf) * o_it->second->currentTransformation ());
	  globalOC_ = M.transform (o_it->first.center ());
          value_type dp = f_it->first.distance (f_it->first.intersection
					  (globalOC_, f_it->first.normal ())),
	    dn = f_it->first.normal ().dot (globalOC_ - f_it->first.center ());
          if (dp < 0) {
	    isInside_ = true;
	    dist = dn * dn;
	  }
          else {
            dist = dp*dp + dn * dn;
	    isInside_ = false;
	  }

          if (dist < minDist) {
            minDist = dist;
            object = o_it;
            floor = f_it;
          }
        }
      }
      relativeTransformation_->joint1 (floor->second);
      relativeTransformation_->joint2 (object->second);
      relativeTransformation_->frame1inJoint1
	(inverse (floor->first.inversePosition ()));
      relativeTransformation_->frame2inJoint2
	(inverse (object->first.inversePosition ()));
    }

    StaticStabilityGravityComplement::StaticStabilityGravityComplement
    (const std::string& name, const std::string& complementName,
     const DevicePtr_t& robot) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (), 3,
			      complementName),
      sibling_ (StaticStabilityGravity::create (name, robot))
    {
    }

    std::pair < StaticStabilityGravityPtr_t,
		StaticStabilityGravityComplementPtr_t >
    StaticStabilityGravityComplement::createPair
    (const std::string& name, const std::string& complementName,
     const DevicePtr_t& robot)
    {
      StaticStabilityGravityComplement* ptr =
	new StaticStabilityGravityComplement (name, complementName, robot);
      StaticStabilityGravityComplementPtr_t shPtr (ptr);
      return std::make_pair (ptr->sibling_, shPtr);
    }

    void StaticStabilityGravityComplement::impl_compute
    (vectorOut_t result, ConfigurationIn_t argument) const
    {
      vector5_t tmp;
      sibling_->impl_compute (tmp, argument);
      result [0] = sibling_->result_ [1];
      result [1] = sibling_->result_ [2];
      result [2] = sibling_->result_ [3];
      if (sibling_->isInside_) {
	result [0] = 0;
	result [1] = 0;
      }
      hppDout (info, "result = " << result.transpose ());
    }

    void StaticStabilityGravityComplement::impl_jacobian
    (matrixOut_t jacobian, ConfigurationIn_t argument) const
    {
      sibling_->computeInternalJacobian (argument);
      if (sibling_->isInside_) {
	jacobian.row (0) = sibling_->jacobian_.row (1);
	jacobian.row (1) = sibling_->jacobian_.row (2);
      } else {
	jacobian.row (0).setZero ();
	jacobian.row (1).setZero ();
      }
      jacobian.row (2) = sibling_->jacobian_.row (3);
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
