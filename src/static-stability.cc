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
#include <limits>
#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>

namespace hpp {
  namespace constraints {
    static void cross (const fcl::Vec3f& v, Eigen::Matrix <double, 3, 3>& m)
    {
      m (0,1) = -v [2]; m (1,0) =  v [2];
      m (0,2) =  v [1]; m (2,0) = -v [1];
      m (1,2) = -v [0]; m (2,1) =  v [0];
    }

    void fclToEigen (const fcl::Vec3f& v, vectorOut_t res)
    {
      res [0] = v[0]; res [1] = v[1]; res [2] = v[2];
    }

    fcl::Vec3f StaticStabilityGravity::gravity (0,0,-1);

    StaticStabilityGravity::StaticStabilityGravity (const DevicePtr_t& robot,
        const JointPtr_t& joint, const fcl::Vec3f& com):
      DifferentiableFunction (robot->configSize (), robot->numberDof (), 3, "StaticStabilityGravity"),
      robot_ (robot), joint_ (joint), com_ (com)
    {
      n_.resize (3);
    }

    StaticStabilityGravityPtr_t StaticStabilityGravity::create (const DevicePtr_t& robot,
        const JointPtr_t& joint, const fcl::Vec3f& com)
    {
      return StaticStabilityGravityPtr_t (new StaticStabilityGravity (robot, joint, com));
    }

    void StaticStabilityGravity::addObjectTriangle (const Triangle& t)
    {
      objectTriangles_.push_back (t);
    }

    void StaticStabilityGravity::addFloorTriangle (const Triangle& t)
    {
      floorTriangles_.push_back (t);
    }

    void StaticStabilityGravity::impl_compute (vectorOut_t result, ConfigurationIn_t argument) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();

      selectTriangles ();

      const Transform3f& M = joint_->currentTransformation ();
      //result [0] = gravity.dot (object_->center () - floor_->center ());
      result [0] = floor_->normal ().dot (M.transform (object_->center ()) - floor_->center ());
      result [1] = floor_->planeXaxis ().dot (M.getRotation () * object_->normal ());
      result [2] = floor_->planeYaxis ().dot (M.getRotation () * object_->normal ());
    }

    void StaticStabilityGravity::impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();

      selectTriangles ();

      const Transform3f& M = joint_->currentTransformation ();
      const JointJacobian_t& Jjoint (joint_->jacobian ());
      cross (M.getRotation () * object_->center (), Rcx_);
      fclToEigen (floor_->normal (), n_);
      jacobian.row (0).leftCols (Jjoint.cols ()) = n_.transpose () *
        (- Rcx_ * Jjoint.bottomRows (3) + Jjoint.topRows (3));
      cross (M.getRotation () * object_->normal (), Rcx_);
      fclToEigen (floor_->planeXaxis (), n_);
      jacobian.row (1).leftCols (Jjoint.cols ()) = - n_.transpose () * Rcx_ * Jjoint.bottomRows (3);
      fclToEigen (floor_->planeYaxis (), n_);
      jacobian.row (2).leftCols (Jjoint.cols ()) = - n_.transpose () * Rcx_ * Jjoint.bottomRows (3);
    }

    void StaticStabilityGravity::selectTriangles () const
    {
      const Transform3f& M = joint_->currentTransformation ();
      // Select the object triangle.
      value_type scalar, maxScalar = - std::numeric_limits <value_type>::infinity();
      for (Triangles::const_iterator it = objectTriangles_.begin ();
          it != objectTriangles_.end (); it++) {
        scalar = gravity.dot (M.getRotation () * it->normal ());
        if (scalar > maxScalar) {
          maxScalar = scalar;
          object_ = it;
        }
      }
      const fcl::Vec3f& c = object_->center ();

      /// Then, select the triangle vertically aligned with this object triangle
      bool hasOneInside = false;
      comGlobalFrame_ = M.transform (com_);
      std::vector <Triangles::const_iterator> selected;
      for (Triangles::const_iterator i = floorTriangles_.begin ();
          i != floorTriangles_.end (); i++) {
        if (i->isInside (comGlobalFrame_, gravity)) {
          selected.push_back (i);
          hasOneInside = true;
        }
      }

      /// Then select the closest triangle along gravity axis.
      assert (hasOneInside);
      maxScalar = - std::numeric_limits <value_type>::infinity();
      for (size_t i = 0; i < selected.size (); i++) {
        scalar = gravity.dot (selected [i]->center () - c);
        if (std::abs (scalar) > maxScalar) {
          maxScalar = std::abs (scalar);
          floor_ = selected [i];
        }
      }
    }
  } // namespace constraints
} // namespace hpp
