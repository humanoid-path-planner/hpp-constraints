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
    static void cross (const fcl::Vec3f& v, eigen::matrix3_t& m)
    {
      m (0,1) = -v [2]; m (1,0) =  v [2];
      m (0,2) =  v [1]; m (2,0) = -v [1];
      m (1,2) = -v [0]; m (2,1) =  v [0];
      m (0,0) = m (1,1) = m (2,2) = 0;
    }

    void fclToEigen (const fcl::Vec3f& v, vectorOut_t res)
    {
      res [0] = v[0]; res [1] = v[1]; res [2] = v[2];
    }

    fcl::Vec3f StaticStabilityGravity::gravity (0,0,-1);

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
      objectTriangles_.push_back (Triangle(t));
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

      const Transform3f& M = joint_->currentTransformation ();
      fcl::Vec3f p = floor_->transformation ().transform (M.transform (object_->center ()));
      result [0] = p[0];
      if (floor_->isInside (floor_->intersection (M.transform (object_->center ()), floor_->normal ()))) {
        result [1] = 0;
        result [2] = 0;
      } else {
        result [1] = p[1];
        result [2] = p[2];
      }

      Rerror_ = object_->transformation ().getRotation () *
        (floor_->transformation ().getRotation () * M.getRotation ()).transpose ();
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

      const Transform3f& M = joint_->currentTransformation ();
      const JointJacobian_t& Jjoint (joint_->jacobian ());
      cross (M.getRotation () * object_->center (), Rcx_);
      eigen::matrix3_t RTeigen;
      do {
        const fcl::Matrix3f& R = floor_->transformation ().getRotation ();
        RTeigen << R (0,0), R (0, 1), R (0, 2),
                   R (1,0), R (1, 1), R (1, 2),
                   R (2,0), R (2, 1), R (2, 2);
      } while (0);
      jacobian_.leftCols (Jjoint.cols ()) = RTeigen *
        (- Rcx_ * Jjoint.bottomRows (3) + Jjoint.topRows (3));
      jacobian.row (0) = jacobian_.row (0);
      if (floor_->isInside (floor_->intersection (M.transform (object_->center ()), floor_->normal ()))) {
        jacobian.row (1).setZero ();
        jacobian.row (2).setZero ();
      } else {
        jacobian.row (1) = jacobian_.row (1);
        jacobian.row (2) = jacobian_.row (2);
      }
      Rerror_ = object_->transformation ().getRotation () * (floor_->transformation ().getRotation () * M.getRotation ()).transpose ();
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
        const fcl::Matrix3f R = object_->transformation ().getRotation () * RT;
        RTeigen << R (0,0), R (0, 1), R (0, 2),
                   R (1,0), R (1, 1), R (1, 2),
                   R (2,0), R (2, 1), R (2, 2);
      } while (0);
      jacobian_.leftCols (Jjoint.cols ()) = -Jlog_ * RTeigen * Jjoint.bottomRows (3);
      jacobian.row (3) = jacobian_.row (1);
      jacobian.row (4) = jacobian_.row (2);
    }

    void StaticStabilityGravity::selectTriangles () const
    {
      const Transform3f& M = joint_->currentTransformation ();
      fcl::Vec3f globalOC_;

      value_type dist, minDist = + std::numeric_limits <value_type>::infinity();
      for (Triangles::const_iterator o_it = objectTriangles_.begin ();
          o_it != objectTriangles_.end (); ++o_it) {
        globalOC_ = M.transform (o_it->center ());
        for (Triangles::const_iterator f_it = floorTriangles_.begin ();
            f_it != floorTriangles_.end (); ++f_it) {
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
  } // namespace constraints
} // namespace hpp
