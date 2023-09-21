// Copyright (c) 2014, LAAS-CNRS
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include "hpp/constraints/convex-shape-contact.hh"

#include <../src/generic-transformation/helper.hh>
#include <hpp/constraints/matrix-view.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/liegroup-element.hh>
#include <limits>

namespace hpp {
namespace constraints {

ConvexShapeContact::ConvexShapeContact(const std::string& name,
                                       DevicePtr_t robot,
                                       const JointAndShapes_t& floorSurfaces,
                                       const JointAndShapes_t& objectSurfaces)
    : DifferentiableFunction(robot->configSize(), robot->numberDof(),
                             LiegroupSpace::Rn(5), name),
      robot_(robot),
      relativeTransformationModel_(robot->numberDof() -
                                   robot->extraConfigSpace().dimension()),
      normalMargin_(0),
      M_(0) {
  relativeTransformationModel_.fullPos = true;
  relativeTransformationModel_.fullOri = true;
  relativeTransformationModel_.rowOri = 3;

  activeParameters_.setConstant(false);
  activeDerivativeParameters_.setConstant(false);
  // Register convex polygons into constraint
  for (JointAndShapes_t::const_iterator it(floorSurfaces.begin());
       it != floorSurfaces.end(); ++it) {
    addFloor(ConvexShape(it->second, it->first));
  }
  for (JointAndShapes_t::const_iterator it(objectSurfaces.begin());
       it != objectSurfaces.end(); ++it) {
    addObject(ConvexShape(it->second, it->first));
  }
  computeRadius();
}

bool ConvexShapeContact::isEqual(const DifferentiableFunction& other) const {
  const ConvexShapeContact& castother =
      dynamic_cast<const ConvexShapeContact&>(other);
  if (!DifferentiableFunction::isEqual(other)) return false;
  if (robot_ != castother.robot_) return false;
  if (objectConvexShapes_.size() != castother.objectConvexShapes_.size())
    return false;
  for (std::size_t i = 0; i < objectConvexShapes_.size(); i++)
    if (floorConvexShapes_[i] != castother.floorConvexShapes_[i]) return false;
  return true;
}

ConvexShapeContactPtr_t ConvexShapeContact::create(
    const std::string& name, DevicePtr_t robot,
    const JointAndShapes_t& floorSurfaces,
    const JointAndShapes_t& objectSurfaces) {
  return ConvexShapeContactPtr_t(
      new ConvexShapeContact(name, robot, floorSurfaces, objectSurfaces));
}

void ConvexShapeContact::addObject(const ConvexShape& t) {
  objectConvexShapes_.push_back(t);

  for (ConvexShapes_t::const_iterator f_it = floorConvexShapes_.begin();
       f_it != floorConvexShapes_.end(); ++f_it) {
    setActiveParameters(robot_, f_it->joint_, t.joint_, activeParameters_,
                        activeDerivativeParameters_);
  }
}

void ConvexShapeContact::addFloor(const ConvexShape& t) {
  ConvexShape tt(t);
  tt.reverse();
  floorConvexShapes_.push_back(tt);

  for (ConvexShapes_t::const_iterator o_it = objectConvexShapes_.begin();
       o_it != objectConvexShapes_.end(); ++o_it) {
    setActiveParameters(robot_, tt.joint_, o_it->joint_, activeParameters_,
                        activeDerivativeParameters_);
  }
}

void ConvexShapeContact::computeRadius() {
  // Compute upper bound of distance between center of polygon and
  // vectices for all floor polygons.
  for (ConvexShapes_t::const_iterator shape(floorConvexShapes_.begin());
       shape != floorConvexShapes_.end(); ++shape) {
    for (std::vector<vector3_t>::const_iterator itv(shape->Pts_.begin());
         itv != shape->Pts_.end(); ++itv) {
      value_type r((*itv - shape->C_).norm());
      if (r > M_) {
        M_ = r;
      }
    }
  }
  M_ += 1;
}

void ConvexShapeContact::setNormalMargin(const value_type& margin) {
  assert(margin >= 0);
  normalMargin_ = margin;
}

std::vector<ConvexShapeContact::ForceData>
ConvexShapeContact::computeContactPoints(ConfigurationIn_t q,
                                         const value_type& normalMargin) const {
  pinocchio::DeviceSync device(robot_);
  device.currentConfiguration(q);
  device.computeForwardKinematics(pinocchio::JOINT_POSITION);

  std::vector<ForceData> forceDatas;
  ForceData forceData;
  ConvexShapeData od, fd;
  for (ConvexShapes_t::const_iterator o_it = objectConvexShapes_.begin();
       o_it != objectConvexShapes_.end(); ++o_it) {
    od.updateToCurrentTransform(*o_it, device.d());
    for (ConvexShapes_t::const_iterator f_it = floorConvexShapes_.begin();
         f_it != floorConvexShapes_.end(); ++f_it) {
      fd.updateToCurrentTransform(*f_it, device.d());
      if (fd.isInside(*f_it, od.center_, fd.normal_)) {
        value_type dn = fd.normal_.dot(od.center_ - fd.center_);
        if (dn < normalMargin) {
          // TODO: compute which points of the object are inside the floor
          // shape.
          forceData.joint = o_it->joint_;
          forceData.points = o_it->Pts_;
          forceData.normal = f_it->N_;
          forceData.supportJoint = f_it->joint_;
          forceDatas.push_back(forceData);
        }
      }
    }
  }
  return forceDatas;
}

void ConvexShapeContact::computeInternalValue(const ConfigurationIn_t& argument,
                                              bool& isInside, ContactType& type,
                                              vector6_t& value,
                                              std::size_t& iobject,
                                              std::size_t& ifloor) const {
  GTDataV<true, true, true, false> data(relativeTransformationModel_, robot_);

  data.device.currentConfiguration(argument);
  data.device.computeForwardKinematics(pinocchio::JOINT_POSITION);

  isInside = selectConvexShapes(data.device.d(), iobject, ifloor);
  const ConvexShape &object(objectConvexShapes_[iobject]),
      floor(floorConvexShapes_[ifloor]);
  type = contactType(object, floor);

  relativeTransformationModel_.joint1 = floor.joint_;
  relativeTransformationModel_.joint2 = object.joint_;
  relativeTransformationModel_.F1inJ1 = floor.positionInJoint();
  relativeTransformationModel_.F2inJ2 = object.positionInJoint();
  relativeTransformationModel_.checkIsIdentity1();
  relativeTransformationModel_.checkIsIdentity2();

  compute<true, true, true, false>::error(data);
  value = data.value;
}

void ConvexShapeContact::impl_compute(LiegroupElementRef result,
                                      ConfigurationIn_t argument) const {
  bool isInside;
  ContactType type;
  vector6_t value;
  std::size_t iobject, ifloor;
  computeInternalValue(argument, isInside, type, value, iobject, ifloor);

  if (isInside) {
    result.vector()[0] = value[0] + normalMargin_;
    result.vector().segment<2>(1).setZero();
  } else {
    result.vector().segment<3>(0) = value.head<3>();
    result.vector()[0] += normalMargin_;
  }
  switch (type) {
    case POINT_ON_PLANE:
      result.vector().segment<2>(3).setZero();
      break;
    case LINE_ON_PLANE:
      // FIXME: only one rotation should be constrained in that case but
      // the relative transformation is not aligned properly. The Y-axis
      // of the reference of "object" should be aligned with the
      // "floor" line axis (Y-axis) projection onto the plane plane.
      // result [3] = 0;
      // result [4] = rt_res_lge[5];
    case PLANE_ON_PLANE:
      result.vector().segment<2>(3) = value.tail<2>();
      break;
  }
  hppDout(info, "result = " << result);
}

void ConvexShapeContact::computeInternalJacobian(
    const ConfigurationIn_t& argument, bool& isInside, ContactType& type,
    matrix_t& jacobian) const {
  static std::vector<bool> mask(6, true);

  GTDataJ<true, true, true, false> data(relativeTransformationModel_, robot_);

  data.device.currentConfiguration(argument);
  data.device.computeForwardKinematics(pinocchio::JOINT_POSITION | pinocchio::JACOBIAN);

  std::size_t ifloor, iobject;
  isInside = selectConvexShapes(data.device.d(), iobject, ifloor);
  const ConvexShape &object(objectConvexShapes_[iobject]),
      floor(floorConvexShapes_[ifloor]);
  type = contactType(object, floor);

  relativeTransformationModel_.joint1 = floor.joint_;
  relativeTransformationModel_.joint2 = object.joint_;
  relativeTransformationModel_.F1inJ1 = floor.positionInJoint();
  relativeTransformationModel_.F2inJ2 = object.positionInJoint();
  relativeTransformationModel_.checkIsIdentity1();
  relativeTransformationModel_.checkIsIdentity2();
  data.cross2.setZero();

  compute<true, true, true, false>::error(data);
  compute<true, true, true, false>::jacobian(data, jacobian, mask);
}

void ConvexShapeContact::impl_jacobian(matrixOut_t jacobian,
                                       ConfigurationIn_t argument) const {
  bool isInside;
  ContactType type;
  matrix_t tmpJac(6, robot_->numberDof());
  computeInternalJacobian(argument, isInside, type, tmpJac);

  if (isInside) {
    jacobian.row(0) = tmpJac.row(0);
    jacobian.row(1).setZero();
    jacobian.row(2).setZero();
  } else {
    jacobian.topRows<3>() = tmpJac.topRows<3>();
  }
  switch (type) {
    case POINT_ON_PLANE:
      jacobian.bottomRows<2>().setZero();
      break;
    case LINE_ON_PLANE:
      // FIXME: See FIXME of impl_compute
      // jacobian.row (3).setZero ();
      // jacobian.row (4) = tmpJac.row (5);
      throw std::logic_error("Contact LINE_ON_PLANE: Unimplement feature");
    case PLANE_ON_PLANE:
      //             Row: 3 4                     Row:  4 5
      jacobian.bottomRows<2>() = tmpJac.bottomRows<2>();
      break;
  }
}

bool ConvexShapeContact::selectConvexShapes(const pinocchio::DeviceData& data,
                                            std::size_t& iobject,
                                            std::size_t& ifloor) const {
  ConvexShapeData od, fd;
  bool isInside = false;  // Initialized only to remove compiler warning.

  value_type dist, minDist = +std::numeric_limits<value_type>::infinity();
  for (std::size_t j = 0; j < floorConvexShapes_.size(); ++j) {
    fd.updateToCurrentTransform(floorConvexShapes_[j], data);

    for (std::size_t i = 0; i < objectConvexShapes_.size(); ++i) {
      od.updateToCurrentTransform(objectConvexShapes_[i], data);
      value_type dp = fd.distance(floorConvexShapes_[j],
                                  fd.intersection(od.center_, fd.normal_)),
                 dn = fd.normal_.dot(od.center_ - fd.center_);
      if (dp < 0)
        dist = dn * dn;
      else
        dist = dp * dp + dn * dn;

      if (dist < minDist) {
        minDist = dist;
        iobject = i;
        ifloor = j;
        isInside = (dp < 0);
      }
    }
  }

  return isInside;
}

ConvexShapeContact::ContactType ConvexShapeContact::contactType(
    const ConvexShape& object, const ConvexShape& floor) const {
  assert(floor.shapeDimension_ > 0 && object.shapeDimension_);
  switch (floor.shapeDimension_) {
    case 1:
      throw std::logic_error("Contact on points is currently unimplemented");
      break;
    case 2:
      throw std::logic_error("Contact on lines is currently unimplemented");
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

ConvexShapeContactComplement::ConvexShapeContactComplement(
    const std::string& name, DevicePtr_t robot,
    const JointAndShapes_t& floorSurfaces,
    const JointAndShapes_t& objectSurfaces)
    : DifferentiableFunction(robot->configSize(), robot->numberDof(),
                             LiegroupSpace::Rn(3),
                             name + std::string("/complement")),
      sibling_(ConvexShapeContact::create(name, robot, floorSurfaces,
                                          objectSurfaces)) {}

bool ConvexShapeContactComplement::isEqual(
    const DifferentiableFunction& other) const {
  const ConvexShapeContactComplement& castother =
      dynamic_cast<const ConvexShapeContactComplement&>(other);
  if (!DifferentiableFunction::isEqual(other)) return false;
  if (sibling_ == castother.sibling_) return true;
  return (*sibling_ == *(castother.sibling_));
}

std::pair<ConvexShapeContactPtr_t, ConvexShapeContactComplementPtr_t>
ConvexShapeContactComplement::createPair(
    const std::string& name, DevicePtr_t robot,
    const JointAndShapes_t& floorSurfaces,
    const JointAndShapes_t& objectSurfaces) {
  ConvexShapeContactComplement* ptr = new ConvexShapeContactComplement(
      name, robot, floorSurfaces, objectSurfaces);
  ConvexShapeContactComplementPtr_t shPtr(ptr);
  return std::make_pair(ptr->sibling_, shPtr);
}

void ConvexShapeContactComplement::computeRelativePoseRightHandSide(
    LiegroupElementConstRef rhs, std::size_t& ifloor, std::size_t& iobject,
    LiegroupElementRef relativePoseRhs) const {
  assert(*(rhs.space()) == *LiegroupSpace::Rn(8));
  value_type M(sibling_->radius());
  vector8_t logRhs(log(rhs));
  vector6_t logRelativePoseRhs;
  logRelativePoseRhs.fill(sqrt(-1));
  // x
  logRelativePoseRhs[0] = logRhs[0];
  // ry, rz
  logRelativePoseRhs.segment<2>(4) = logRhs.segment<2>(3);
  // rx
  logRelativePoseRhs[3] = logRhs[7];
  value_type Y(logRhs[5]);
  value_type Z(logRhs[6]);
  assert(floor(Y / (2 * M) + .5) >= 0);
  ifloor = (std::size_t)(floor(Y / (2 * M) + .5));
  assert(floor(Z / (2 * M) + .5) >= 0);
  iobject = (std::size_t)(floor(Z / (2 * M) + .5));
  if ((logRhs[1] == 0) && (logRhs[2] == 0)) {
    // inside
    logRelativePoseRhs[1] = Y - (value_type)(2 * ifloor) * M;
    logRelativePoseRhs[2] = Z - (value_type)(2 * iobject) * M;
  } else {
    // outside
    logRelativePoseRhs.segment<2>(1) = logRhs.segment<2>(1);
  }
  assert(!logRelativePoseRhs.hasNaN());
  relativePoseRhs = LiegroupSpace::R3xSO3()->exp(logRelativePoseRhs);
  assert(*(relativePoseRhs.space()) == *LiegroupSpace::R3xSO3());
}

void ConvexShapeContactComplement::impl_compute(
    LiegroupElementRef result, ConfigurationIn_t argument) const {
  value_type M(sibling_->radius());
  bool isInside;
  ConvexShapeContact::ContactType type;
  vector6_t value;
  std::size_t iobject, ifloor;
  sibling_->computeInternalValue(argument, isInside, type, value, iobject,
                                 ifloor);
  vector3_t offset;
  offset << (value_type)(2 * ifloor) * M, (value_type)(2 * iobject) * M, 0;
  result.vector()[2] = value[3];
  if (isInside)
    result.vector().head<2>() = value.segment<2>(1);
  else
    result.vector().head<2>().setZero();
  result.vector() += offset;
  hppDout(info, "result = " << result);
}

void ConvexShapeContactComplement::impl_jacobian(
    matrixOut_t jacobian, ConfigurationIn_t argument) const {
  bool isInside;
  ConvexShapeContact::ContactType type;
  matrix_t tmpJac(6, sibling_->robot_->numberDof());
  sibling_->computeInternalJacobian(argument, isInside, type, tmpJac);

  if (isInside)
    jacobian.topRows<2>() = tmpJac.middleRows<2>(1);
  else
    jacobian.topRows<2>().setZero();
  jacobian.row(2) = tmpJac.row(3);
}

std::ostream& ConvexShapeContact::print(std::ostream& o) const {
  o << "ConvexShapeContact: " << name() << ", active dof "
    << pretty_print(BlockIndex::fromLogicalExpression(activeParameters_))
    << incindent;

  o << iendl << "Object shapes:" << incindent << iendl;
  for (ConvexShapes_t::const_iterator o_it = objectConvexShapes_.begin();
       o_it != objectConvexShapes_.end(); ++o_it) {
    if (o_it->joint_)
      o << "object on " << o_it->joint_->name() << iendl;
    else
      o << "object on universe" << iendl;
    o << "position in joint:" << iendl;
    o << incindent << o_it->positionInJoint();
    o << decindent << iendl;
  }
  for (ConvexShapes_t::const_iterator fl_it = floorConvexShapes_.begin();
       fl_it != floorConvexShapes_.end(); ++fl_it) {
    if (fl_it->joint_)
      o << "floor on " << fl_it->joint_->name() << iendl;
    else
      o << "floor on universe" << iendl;
    o << "position in joint:" << iendl;
    o << incindent << fl_it->positionInJoint();
    o << decindent << iendl;
  }
  return o << decindent;
}

ConvexShapeContactHoldPtr_t ConvexShapeContactHold::create(
    const std::string& name, DevicePtr_t robot,
    const JointAndShapes_t& floorSurfaces,
    const JointAndShapes_t& objectSurfaces) {
  ConvexShapeContactHold* ptr(
      new ConvexShapeContactHold(name, robot, floorSurfaces, objectSurfaces));
  return ConvexShapeContactHoldPtr_t(ptr);
}

ConvexShapeContactHold::ConvexShapeContactHold(
    const std::string& name, DevicePtr_t robot,
    const JointAndShapes_t& floorSurfaces,
    const JointAndShapes_t& objectSurfaces)
    : DifferentiableFunction(robot->configSize(), robot->numberDof(),
                             LiegroupSpace::Rn(8), name) {
  std::pair<ConvexShapeContactPtr_t, ConvexShapeContactComplementPtr_t> pair(
      ConvexShapeContactComplement::createPair(name, robot, floorSurfaces,
                                               objectSurfaces));
  constraint_ = pair.first;
  complement_ = pair.second;
}

void ConvexShapeContactHold::impl_compute(LiegroupElementRef result,
                                          vectorIn_t argument) const {
  LiegroupElementRef tmp1(result.vector().head<5>(), LiegroupSpace::Rn(5));
  LiegroupElementRef tmp2(result.vector().tail<3>(), LiegroupSpace::Rn(3));
  constraint_->impl_compute(tmp1, argument);
  complement_->impl_compute(tmp2, argument);
}

void ConvexShapeContactHold::impl_jacobian(matrixOut_t jacobian,
                                           vectorIn_t arg) const {
  constraint_->impl_jacobian(jacobian.topRows<5>(), arg);
  complement_->impl_jacobian(jacobian.bottomRows<3>(), arg);
}

bool ConvexShapeContactHold::isEqual(
    const DifferentiableFunction& other) const {
  const ConvexShapeContactHold& castother =
      dynamic_cast<const ConvexShapeContactHold&>(other);
  if (!DifferentiableFunction::isEqual(other)) return false;
  if (constraint_ != castother.constraint_) return false;
  if (complement_ != castother.complement_) return false;
  return true;
}

}  // namespace constraints
}  // namespace hpp
