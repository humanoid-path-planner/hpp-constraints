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

#include "hpp/constraints/convex-shape-contact.hh"

#include <limits>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/liegroup-element.hh>

#include <hpp/constraints/matrix-view.hh>

namespace hpp {
  namespace constraints {

    ConvexShapeContact::ConvexShapeContact
    (const std::string& name, const DevicePtr_t& robot) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
                              LiegroupSpace::Rn (5), name), robot_ (robot),
      relativeTransformation_ (name, robot, std::vector<bool>(6, true)),
      normalMargin_ (0), result_ (LiegroupSpace::Rn (6))
    {
      relativeTransformation_.joint1(robot->rootJoint());
      relativeTransformation_.joint2(robot->rootJoint());
      jacobian_.resize (6, robot->numberDof ());
      activeParameters_.setConstant(false);
      activeDerivativeParameters_.setConstant(false);
    }

    ConvexShapeContactPtr_t ConvexShapeContact::create (
        const std::string& name,
        const DevicePtr_t& robot)
    {
      return ConvexShapeContactPtr_t (new ConvexShapeContact
					  (name, robot));
    }

    ConvexShapeContactPtr_t ConvexShapeContact::create
    (const DevicePtr_t& robot)
    {
      return create ("ConvexShapeContact", robot);
    }

    void ConvexShapeContact::addObjectTriangle (const fcl::TriangleP& t,
						    const JointPtr_t& joint)
    {
      addObject (ConvexShape (t, joint));
    }

    void ConvexShapeContact::addFloorTriangle (const fcl::TriangleP& t,
						   const JointPtr_t& joint)
    {
      addFloor (ConvexShape (t, joint));
    }

    void ConvexShapeContact::addObject (const ConvexShape& t)
    {
      objectConvexShapes_.push_back (t);

      relativeTransformation_.joint2 (t.joint_);
      for (ConvexShapes_t::const_iterator f_it = floorConvexShapes_.begin ();
          f_it != floorConvexShapes_.end (); ++f_it) {
        relativeTransformation_.joint1 (f_it->joint_);
        activeParameters_ = activeParameters_
          || relativeTransformation_.activeParameters();
        activeDerivativeParameters_ = activeDerivativeParameters_
          || relativeTransformation_.activeDerivativeParameters();
      }
    }

    void ConvexShapeContact::addFloor (const ConvexShape& t)
    {
      ConvexShape tt (t); tt.reverse ();
      floorConvexShapes_.push_back (tt);

      relativeTransformation_.joint1 (tt.joint_);
      for (ConvexShapes_t::const_iterator o_it = objectConvexShapes_.begin ();
          o_it != objectConvexShapes_.end (); ++o_it) {
        relativeTransformation_.joint2 (o_it->joint_);
        activeParameters_ = activeParameters_
          || relativeTransformation_.activeParameters();
        activeDerivativeParameters_ = activeDerivativeParameters_
          || relativeTransformation_.activeDerivativeParameters();
      }
    }

    void ConvexShapeContact::setNormalMargin (const value_type& margin)
    {
      assert (margin >= 0);
      normalMargin_ = margin;
    }

    std::vector <ConvexShapeContact::ForceData>
      ConvexShapeContact::computeContactPoints (
          const value_type& normalMargin) const
    {
      std::vector <ForceData> fds;
      ForceData fd;
      for (ConvexShapes_t::const_iterator o_it = objectConvexShapes_.begin ();
          o_it != objectConvexShapes_.end (); ++o_it) {
        // o_it->updateToCurrentTransform ();
        const vector3_t& globalOC_ = o_it->center ();
        for (ConvexShapes_t::const_iterator f_it = floorConvexShapes_.begin ();
            f_it != floorConvexShapes_.end (); ++f_it) {
          // f_it->updateToCurrentTransform ();
          if (f_it->isInside (globalOC_, f_it->normal ())) {
            value_type dn = f_it->normal ().dot (globalOC_ - f_it->center ());
            if (dn < normalMargin) {
              // TODO: compute which points of the object are inside the floor shape.
              fd.joint = o_it->joint_;
              fd.points = o_it->Pts_;
              fd.normal = f_it->N_;
              fd.supportJoint = f_it->joint_;
              fds.push_back (fd);
            }
          }
        }
      }
      return fds;
    }

    void ConvexShapeContact::impl_compute (LiegroupElement& result,
                                           ConfigurationIn_t argument) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();

      selectConvexShapes ();
      relativeTransformation_.value (result_, argument);
      if (isInside_) {
        result.vector () [0] = result_.vector () [0] + normalMargin_;
        result.vector ().segment <2> (1).setZero ();
      } else {
        result.vector ().segment <3> (0) = result_.vector ().segment <3> (0);
        result.vector () [0] += normalMargin_;
      }
      switch (contactType_) {
        case POINT_ON_PLANE:
          result.vector ().segment <2> (3).setZero ();
          break;
        case LINE_ON_PLANE:
          // FIXME: only one rotation should be constrained in that case but
          // the relative transformation is not aligned properly. The Y-axis
          // of the reference of "object" should be aligned with the
          // "floor" line axis (Y-axis) projection onto the plane plane.
          // result [3] = 0;
          // result [4] = result_[5];
        case PLANE_ON_PLANE:
          result.vector ().segment<2> (3) = result_.vector ().segment<2> (4);
          break;
      }
      hppDout (info, "result = " << result);
    }

    void ConvexShapeContact::computeInternalJacobian
    (ConfigurationIn_t argument) const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      selectConvexShapes ();
      relativeTransformation_.jacobian (jacobian_, argument);
    }

    void ConvexShapeContact::impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const
    {
      computeInternalJacobian (argument);
      if (isInside_) {
	jacobian.row (0) = jacobian_.row (0);
	jacobian.row (1).setZero ();
	jacobian.row (2).setZero ();
      } else {
        jacobian.topRows<3> () = jacobian_.topRows <3> ();
      }
      switch (contactType_) {
        case POINT_ON_PLANE:
          jacobian.bottomRows<2> ().setZero ();
          break;
        case LINE_ON_PLANE:
          // FIXME: See FIXME of impl_compute
          // jacobian.row (3).setZero ();
          // jacobian.row (4) = jacobian_.row (5);
          throw std::logic_error ("Contact LINE_ON_PLANE: Unimplement feature");
        case PLANE_ON_PLANE:
          //             Row: 3 4                     Row:  4 5
          jacobian.bottomRows<2> () = jacobian_.bottomRows <2> ();
          break;
      }
    }

    void ConvexShapeContact::selectConvexShapes () const
    {
      ConvexShapes_t::const_iterator object;
      ConvexShapes_t::const_iterator floor;

      value_type dist, minDist = + std::numeric_limits <value_type>::infinity();
      for (ConvexShapes_t::const_iterator o_it = objectConvexShapes_.begin ();
          o_it != objectConvexShapes_.end (); ++o_it) {
        o_it->updateToCurrentTransform ();
        const vector3_t& globalOC_ = o_it->center ();
        for (ConvexShapes_t::const_iterator f_it = floorConvexShapes_.begin ();
            f_it != floorConvexShapes_.end (); ++f_it) {
          f_it->updateToCurrentTransform ();
          value_type dp = f_it->distance (f_it->intersection (globalOC_, f_it->normal ())),
                     dn = f_it->normal ().dot (globalOC_ - f_it->center ());
          if (dp < 0) dist = dn * dn;
          else        dist = dp*dp + dn * dn;

          if (dist < minDist) {
            minDist = dist;
            object = o_it;
            floor = f_it;
	    isInside_ = (dp < 0);
          }
        }
      }
      contactType_ = contactType (*object, *floor);
      relativeTransformation_.joint1 (floor->joint_);
      relativeTransformation_.joint2 (object->joint_);
      relativeTransformation_.frame1InJoint1 (floor->positionInJoint ());
      relativeTransformation_.frame2InJoint2 (object->positionInJoint ());
    }

    ConvexShapeContact::ContactType ConvexShapeContact::contactType (
        const ConvexShape& object, const ConvexShape& floor) const
    {
      assert (floor.shapeDimension_ > 0 && object.shapeDimension_);
      switch (floor.shapeDimension_) {
        case 1:
          throw std::logic_error
            ("Contact on points is currently unimplemented");
          break;
        case 2:
          throw std::logic_error
            ("Contact on lines is currently unimplemented");
          break;
        default:
          switch (object.shapeDimension_) {
            case 1:
              return POINT_ON_PLANE;
              break;
            case 2:
              return LINE_ON_PLANE;
              break;
            default:
              return PLANE_ON_PLANE;
              break;
          }
          break;
      }
    }

    ConvexShapeContactComplement::ConvexShapeContactComplement
    (const std::string& name, const std::string& complementName,
     const DevicePtr_t& robot) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (), 3,
			      complementName),
      sibling_ (ConvexShapeContact::create (name, robot))
    {
    }

    std::pair < ConvexShapeContactPtr_t,
		ConvexShapeContactComplementPtr_t >
    ConvexShapeContactComplement::createPair
    (const std::string& name, const std::string& complementName,
     const DevicePtr_t& robot)
    {
      ConvexShapeContactComplement* ptr =
	new ConvexShapeContactComplement (name, complementName, robot);
      ConvexShapeContactComplementPtr_t shPtr (ptr);
      return std::make_pair (ptr->sibling_, shPtr);
    }

    void ConvexShapeContactComplement::impl_compute
    (LiegroupElement& result, ConfigurationIn_t argument) const
    {
      LiegroupElement tmp (LiegroupSpace::Rn (5));
      sibling_->impl_compute (tmp, argument);
      result.vector () [2] = sibling_->result_.vector () [3];
      if (sibling_->isInside_) {
	result.vector () [0] = sibling_->result_.vector () [1];
	result.vector () [1] = sibling_->result_.vector () [2];
      } else {
	result.vector () [0] = 0;
	result.vector () [1] = 0;
      }
      hppDout (info, "result = " << result);
    }

    void ConvexShapeContactComplement::impl_jacobian
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

    std::ostream& ConvexShapeContact::print (std::ostream& o) const
    {
      o << "ConvexShapeContact: " << name () << ", active dof "
        << pretty_print (BlockIndex::fromLogicalExpression (activeParameters_)) << incindent;

      o << iendl << "Object shapes:" << incindent;
      for (ConvexShapes_t::const_iterator o_it = objectConvexShapes_.begin ();
          o_it != objectConvexShapes_.end (); ++o_it) {
        if (o_it->joint_)
          o << "object on " << o_it->joint_->name() << iendl;
        else
          o << "object on universe" << iendl;
      }
      for (ConvexShapes_t::const_iterator fl_it = floorConvexShapes_.begin ();
          fl_it != floorConvexShapes_.end (); ++fl_it) {
        if (fl_it->joint_)
          o << "floor on " << fl_it->joint_->name() << iendl;
        else
          o << "floor on universe" << iendl;
      }
      return o << decindent;
    }
  } // namespace constraints
} // namespace hpp
