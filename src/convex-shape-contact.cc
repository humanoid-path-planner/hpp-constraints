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

#include <../src/generic-transformation/helper.hh>

namespace hpp {
  namespace constraints {

    ConvexShapeContact::ConvexShapeContact
    (const std::string& name, const DevicePtr_t& robot) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
                              LiegroupSpace::Rn (5), name), robot_ (robot),
      relativeTransformationModel_ (robot->numberDof()-robot->extraConfigSpace().dimension()),
      normalMargin_ (0)
    {
      relativeTransformationModel_.fullPos = true;
      relativeTransformationModel_.fullOri = true;
      relativeTransformationModel_.rowOri = 3;

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

      for (ConvexShapes_t::const_iterator f_it = floorConvexShapes_.begin ();
          f_it != floorConvexShapes_.end (); ++f_it) {
        setActiveParameters (robot_, f_it->joint_, t.joint_,
            activeParameters_, activeDerivativeParameters_);
      }
    }

    void ConvexShapeContact::addFloor (const ConvexShape& t)
    {
      ConvexShape tt (t); tt.reverse ();
      floorConvexShapes_.push_back (tt);

      for (ConvexShapes_t::const_iterator o_it = objectConvexShapes_.begin ();
          o_it != objectConvexShapes_.end (); ++o_it) {
        setActiveParameters (robot_, tt.joint_, o_it->joint_,
            activeParameters_, activeDerivativeParameters_);
      }
    }

    void ConvexShapeContact::setNormalMargin (const value_type& margin)
    {
      assert (margin >= 0);
      normalMargin_ = margin;
    }

    std::vector <ConvexShapeContact::ForceData>
      ConvexShapeContact::computeContactPoints (ConfigurationIn_t q,
          const value_type& normalMargin) const
    {
      pinocchio::DeviceSync device (robot_);
      device.currentConfiguration (q);
      device.computeForwardKinematics ();

      std::vector <ForceData> forceDatas;
      ForceData forceData;
      ConvexShapeData od, fd;
      for (ConvexShapes_t::const_iterator o_it = objectConvexShapes_.begin ();
          o_it != objectConvexShapes_.end (); ++o_it) {
        od.updateToCurrentTransform (*o_it, device.d());
        for (ConvexShapes_t::const_iterator f_it = floorConvexShapes_.begin ();
            f_it != floorConvexShapes_.end (); ++f_it) {
          fd.updateToCurrentTransform (*f_it, device.d());
          if (fd.isInside (*f_it, od.center_, fd.normal_)) {
            value_type dn = fd.normal_.dot (od.center_ - fd.center_);
            if (dn < normalMargin) {
              // TODO: compute which points of the object are inside the floor shape.
              forceData.joint = o_it->joint_;
              forceData.points = o_it->Pts_;
              forceData.normal = f_it->N_;
              forceData.supportJoint = f_it->joint_;
              forceDatas.push_back (forceData);
            }
          }
        }
      }
      return forceDatas;
    }

    void ConvexShapeContact::computeInternalValue (const ConfigurationIn_t& argument,
        bool& isInside, ContactType& type, vector6_t& value) const
    {
      GTDataV<true, true, true, false> data (relativeTransformationModel_, robot_);

      data.device.currentConfiguration (argument);
      data.device.computeForwardKinematics ();

      ConvexShapes_t::const_iterator object, floor;
      isInside = selectConvexShapes (data.device.d(), object, floor);
      type = contactType (*object, *floor);

      relativeTransformationModel_.joint1 = floor->joint_;
      relativeTransformationModel_.joint2 = object->joint_;
      relativeTransformationModel_.F1inJ1 = floor->positionInJoint ();
      relativeTransformationModel_.F2inJ2 = object->positionInJoint ();
      relativeTransformationModel_.checkIsIdentity1();
      relativeTransformationModel_.checkIsIdentity2();

      compute<true, true, true, false>::error (data);
      value = data.value;
    }

    void ConvexShapeContact::impl_compute (LiegroupElementRef result,
                                           ConfigurationIn_t argument) const
    {
      bool isInside;
      ContactType type;
      vector6_t value;
      computeInternalValue (argument, isInside, type, value);

      if (isInside) {
        result.vector () [0] = value [0] + normalMargin_;
        result.vector ().segment <2> (1).setZero ();
      } else {
        result.vector ().segment <3> (0) = value.head <3> ();
        result.vector () [0] += normalMargin_;
      }
      switch (type) {
        case POINT_ON_PLANE:
          result.vector ().segment <2> (3).setZero ();
          break;
        case LINE_ON_PLANE:
          // FIXME: only one rotation should be constrained in that case but
          // the relative transformation is not aligned properly. The Y-axis
          // of the reference of "object" should be aligned with the
          // "floor" line axis (Y-axis) projection onto the plane plane.
          // result [3] = 0;
          // result [4] = rt_res_lge[5];
        case PLANE_ON_PLANE:
          result.vector ().segment<2> (3) = value.tail<2> ();
          break;
      }
      hppDout (info, "result = " << result);
    }

    void ConvexShapeContact::computeInternalJacobian
    (const ConfigurationIn_t& argument,
     bool& isInside, ContactType& type, matrix_t& jacobian) const
    {
      static std::vector<bool> mask (6, true);

      GTDataJ<true, true, true, false> data (relativeTransformationModel_, robot_);

      data.device.currentConfiguration (argument);
      data.device.computeForwardKinematics ();

      ConvexShapes_t::const_iterator object, floor;
      isInside = selectConvexShapes (data.device.d(), object, floor);
      type = contactType (*object, *floor);

      relativeTransformationModel_.joint1 = floor->joint_;
      relativeTransformationModel_.joint2 = object->joint_;
      relativeTransformationModel_.F1inJ1 = floor->positionInJoint ();
      relativeTransformationModel_.F2inJ2 = object->positionInJoint ();
      relativeTransformationModel_.checkIsIdentity1();
      relativeTransformationModel_.checkIsIdentity2();
      data.cross2.setZero();

      compute<true, true, true, false>::error (data);
      compute<true, true, true, false>::jacobian (data, jacobian, mask);
    }

    void ConvexShapeContact::impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const
    {
      bool isInside;
      ContactType type;
      matrix_t tmpJac (6, robot_->numberDof());
      computeInternalJacobian (argument, isInside, type, tmpJac);

      if (isInside) {
        jacobian.row (0) = tmpJac.row (0);
        jacobian.row (1).setZero ();
        jacobian.row (2).setZero ();
      } else {
        jacobian.topRows<3> () = tmpJac.topRows <3> ();
      }
      switch (type) {
        case POINT_ON_PLANE:
          jacobian.bottomRows<2> ().setZero ();
          break;
        case LINE_ON_PLANE:
          // FIXME: See FIXME of impl_compute
          // jacobian.row (3).setZero ();
          // jacobian.row (4) = tmpJac.row (5);
          throw std::logic_error ("Contact LINE_ON_PLANE: Unimplement feature");
        case PLANE_ON_PLANE:
          //             Row: 3 4                     Row:  4 5
          jacobian.bottomRows<2> () = tmpJac.bottomRows <2> ();
          break;
      }
    }

    bool ConvexShapeContact::selectConvexShapes (const pinocchio::DeviceData& data,
        ConvexShapes_t::const_iterator& object,
        ConvexShapes_t::const_iterator& floor) const
    {
      ConvexShapeData od, fd;
      bool isInside = false; // Initialized only to remove compiler warning.

      value_type dist, minDist = + std::numeric_limits <value_type>::infinity();
      for (ConvexShapes_t::const_iterator o_it = objectConvexShapes_.begin ();
          o_it != objectConvexShapes_.end (); ++o_it) {
        od.updateToCurrentTransform (*o_it, data);

        for (ConvexShapes_t::const_iterator f_it = floorConvexShapes_.begin ();
            f_it != floorConvexShapes_.end (); ++f_it) {
          fd.updateToCurrentTransform (*f_it, data);
          value_type dp = fd.distance (*f_it, fd.intersection (od.center_, fd.normal_)),
                     dn = fd.normal_.dot (od.center_ - fd.center_);
          if (dp < 0) dist = dn * dn;
          else        dist = dp*dp + dn * dn;

          if (dist < minDist) {
            minDist = dist;
            object = o_it;
            floor = f_it;
            isInside = (dp < 0);
          }
        }
      }

      return isInside;
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
    (LiegroupElementRef result, ConfigurationIn_t argument) const
    {
      bool isInside;
      ConvexShapeContact::ContactType type;
      vector6_t value;
      sibling_->computeInternalValue (argument, isInside, type, value);

      result.vector () [2] = value [3];
      if (isInside) result.vector ().head<2>() = value.segment<2>(1);
      else          result.vector ().head<2>().setZero();
      hppDout (info, "result = " << result);
    }

    void ConvexShapeContactComplement::impl_jacobian
    (matrixOut_t jacobian, ConfigurationIn_t argument) const
    {
      bool isInside;
      ConvexShapeContact::ContactType type;
      matrix_t tmpJac (6, sibling_->robot_->numberDof());
      sibling_->computeInternalJacobian (argument, isInside, type, tmpJac);

      if (isInside)
        jacobian.topRows<2>() = tmpJac.middleRows<2>(1);
      else
        jacobian.topRows<2>().setZero ();
      jacobian.row (2) = tmpJac.row (3);
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
