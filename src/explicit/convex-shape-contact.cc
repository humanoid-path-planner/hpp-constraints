// Copyright (c) 2020, Airbus SAS and CNRS
// Authors: Florent Lamiraux
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

#define HPP_DEBUG
#include <boost/assign/list_of.hpp>
#include <hpp/constraints/convex-shape-contact.hh>
#include <hpp/constraints/explicit/convex-shape-contact.hh>
#include <hpp/constraints/explicit/relative-pose.hh>
#include <hpp/constraints/affine-function.hh>
#include <hpp/pinocchio/liegroup-space.hh>
#include "../src/explicit/input-configurations.hh"

namespace hpp {
  namespace constraints {
    namespace explicit_ {
      using boost::assign::list_of;

      ConvexShapeContactPtr_t ConvexShapeContact::create
      (const std::string& name, DevicePtr_t robot,
       const JointAndShapes_t& floorSurfaces,
       const JointAndShapes_t& objectSurfaces,
       const value_type& margin)
      {
        ConvexShapeContact* ptr(new ConvexShapeContact
                                (name, robot, floorSurfaces, objectSurfaces,
                                 margin));
        ConvexShapeContactPtr_t shPtr(ptr);
        ptr->init(shPtr);
        return shPtr;
      }

      ConvexShapeContact::Constraints_t
      ConvexShapeContact::createConstraintAndComplement
      (const std::string& name, DevicePtr_t robot,
       const JointAndShapes_t& floorSurfaces,
       const JointAndShapes_t& objectSurfaces,
       const value_type& margin)
      {
        Constraints_t result;
        std::pair < hpp::constraints::ConvexShapeContactPtr_t,
                    hpp::constraints::ConvexShapeContactComplementPtr_t >
          functions(ConvexShapeContactComplement::createPair
                    (name, robot, floorSurfaces, objectSurfaces));
        functions.first->setNormalMargin(margin);
        // Contact constraint (= 0)
        std::get<0>(result) = Implicit::create(functions.first, 5*EqualToZero);
        // Contact constraint complement (= rhs)
        std::get<1>(result) = Implicit::create
          (functions.second, ComparisonTypes_t(functions.second->outputSize(),
            constraints::Equality));
        std::get<2>(result) = create(name + "/hold", robot, floorSurfaces,
                                     objectSurfaces, margin);
        return result;
      }

      ConvexShapeContactPtr_t ConvexShapeContact::createCopy
      (const ConvexShapeContactPtr_t& other)
      {
        ConvexShapeContact* ptr (new ConvexShapeContact (*other));
        ConvexShapeContactPtr_t shPtr (ptr);
        ConvexShapeContactWkPtr_t wkPtr (shPtr);
        ptr->init (wkPtr);
        return shPtr;
      }

      ImplicitPtr_t ConvexShapeContact::copy () const
      {
        return createCopy (weak_.lock ());
      }

       void ConvexShapeContact::outputValue
       (LiegroupElementRef result, vectorIn_t qin, vectorIn_t rhs) const
       {
        assert(HPP_DYNAMIC_PTR_CAST(ConvexShapeContactHold, functionPtr()));
        ConvexShapeContactHoldPtr_t f(HPP_STATIC_PTR_CAST
                                      (ConvexShapeContactHold, functionPtr()));
        std::size_t ifloor, iobject;
        vector6_t relativePoseRhs;
        f->complement()->computeRelativePoseRightHandSide
          (rhs, ifloor, iobject, relativePoseRhs);
        // Extract input configuration of relative pose from qin
        Eigen::RowBlockIndices inputIndices(inputConf());
        vector_t q(f->inputSize()); q.fill(sqrt(-1));
        inputIndices.lview(q) = qin;
        RelativePosePtr_t relativePose(pose_[ifloor*nFloor_ + iobject]);
        Eigen::RowBlockIndices relPosInputIndices(relativePose->inputConf());
        vector_t qinRelPose = relPosInputIndices.rview(q);
        assert(!qinRelPose.hasNaN());
        relativePose->outputValue(result, qinRelPose, relativePoseRhs);
       }

      void ConvexShapeContact::jacobianOutputValue
      (vectorIn_t qin, LiegroupElementConstRef, vectorIn_t rhs,
       matrixOut_t jacobian) const
      {
        assert(HPP_DYNAMIC_PTR_CAST(ConvexShapeContactHold, functionPtr()));
        ConvexShapeContactHoldPtr_t f(HPP_STATIC_PTR_CAST
                                      (ConvexShapeContactHold, functionPtr()));
        std::size_t ifloor, iobject;
        vector6_t relativePoseRhs;
        f->complement()->computeRelativePoseRightHandSide
          (rhs, ifloor, iobject, relativePoseRhs);
        // Extract input configuration of relative pose from qin
        Eigen::RowBlockIndices inputIndices(inputConf());
        vector_t q(f->inputSize()); q.fill(sqrt(-1));
        inputIndices.lview(q) = qin;
        RelativePosePtr_t relativePose(pose_[ifloor*nFloor_ + iobject]);
        Eigen::RowBlockIndices relPosInputIndices(relativePose->inputConf());
        vector_t qinRelPose = relPosInputIndices.rview(q);
        assert(!qinRelPose.hasNaN());
        LiegroupElement outputRelPose(LiegroupSpace::SE3());
        relativePose->outputValue(outputRelPose, qinRelPose, relativePoseRhs);
        relativePose->jacobianOutputValue(qinRelPose, outputRelPose,
                                          relativePoseRhs, jacobian);
      }

      ConvexShapeContact::ConvexShapeContact
      (const std::string& name, DevicePtr_t robot,
       const JointAndShapes_t& floorSurfaces,
       const JointAndShapes_t& objectSurfaces,
       const value_type& margin) :
        Explicit(ConvexShapeContactHold::create
                 (name, robot, floorSurfaces, objectSurfaces),
                 ConstantFunction::create
                 (pinocchio::LiegroupSpace::SE3()->neutral(),
                  contact::inputSize(robot, floorSurfaces, objectSurfaces),
                  contact::inputDerivSize(robot,floorSurfaces,objectSurfaces),
                  name),
                 contact::inputConfVariables
                 (robot, floorSurfaces, objectSurfaces),
                 relativePose::jointConfInterval(objectSurfaces.front().first),
                 contact::inputVelocityVariables
                 (robot, floorSurfaces, objectSurfaces),
                 relativePose::jointVelInterval(objectSurfaces.front().first),
                 (5 * EqualToZero << 3 * Equality))
      {
        assert(HPP_DYNAMIC_PTR_CAST(ConvexShapeContactHold, functionPtr()));
        ConvexShapeContactHoldPtr_t f
          (HPP_STATIC_PTR_CAST(ConvexShapeContactHold, functionPtr()));
        f->contactConstraint()->setNormalMargin(margin);
        const ConvexShapes_t& fs
          (f->contactConstraint()->floorContactSurfaces());
        const ConvexShapes_t& os
          (f->contactConstraint()->objectContactSurfaces());
        nFloor_ = fs.size();
        // Compute explicit relative poses
        for (std::size_t j=0; j<fs.size(); ++j)
        {
          // move floor surface along x to take into account margin.
          Transform3f posInJoint(fs[j].positionInJoint());
          hppDout(info, "posInJoint" << posInJoint);
          posInJoint.translation() -= margin * posInJoint.rotation().col(0);
          hppDout(info, "posInJoint" << posInJoint);
          for (std::size_t i=0; i<os.size(); ++i)
          {
            // Create explicit relative pose for each combination
            // (floor surface, object surface)
            pose_.push_back(RelativePose::create
                            ("",robot, fs[j].joint_, os[i].joint_,
                             posInJoint, os[i].positionInJoint(),
                             6 * Equality));
          }
        }
      }
      void ConvexShapeContact::init (ConvexShapeContactWkPtr_t weak)
      {
        Explicit::init (weak);
        weak_ = weak;
      }
    } // namespace explicit_
  } // namespace constraints
} // namespace hpp
