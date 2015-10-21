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

#include <hpp/constraints/contact-model.hh>

#include <hpp/constraints/transformation.hh>
#include <hpp/constraints/relative-transformation.hh>

#define UNIMPLEMENTED throw std::logic_error ("Unimplement feature")

namespace hpp {
  namespace constraints {
    AbstractSurface::Point_t AbstractSurface::pointToGlobalFrame (const Point_t& p) const {
      if (joint_ == NULL) return p;
      else return joint_->currentTransformation ().transform (p);
    }

    AbstractSurface::Point_t AbstractSurface::vectorToGlobalFrame (const Point_t& p) const {
      if (joint_ == NULL) return p;
      else return joint_->currentTransformation ().getRotation () * p;
    }

    AbstractSurface::Point_t AbstractSurface::pointToLocalFrame (const Point_t& p) const {
      if (joint_ == NULL) return p;
      else {
        fcl::Transform3f T = joint_->currentTransformation();
        return T.inverse().transform (p);
      }
    }

    AbstractSurface::Point_t AbstractSurface::vectorToLocalFrame (const Point_t& p) const {
      if (joint_ == NULL) return p;
      else return joint_->currentTransformation ().getRotation ().transposeTimes (p);
    }

    AbstractSurface::ForceDatas_t PointSurface::forcePoints (const AbstractSurface* other) const {
      ForceDatas_t fds (1);
      fds[0].j = joint_;
      fds[0].p.push_back (p_);
      fds[0].n = AbstractSurface::normal (other);
      return fds;
    }

    AbstractSurface::ForceDatas_t LineSurface::forcePoints (const AbstractSurface* other) const {
      ForceDatas_t fds (1);
      fds[0].j = joint_;
      fds[0].p.push_back (p0_);
      fds[0].p.push_back (p1_);
      fds[0].n = AbstractSurface::normal (other);
      return fds;
    }

    AbstractSurface::ForceDatas_t FlatSurface::forcePoints (const AbstractSurface* other) const {
      ForceDatas_t fds (ch_.Pts_.size());
      fds[0].j = joint_;
      fds[0].p = ch_.Pts_;
      fds[0].n = AbstractSurface::normal (other);
      return fds;
    }

    AbstractSurface::Point_t AbstractSurface::normal (const AbstractSurface* other) const {
      const PointSurface* ps = dynamic_cast <const PointSurface*> (other);
      if (ps != NULL) return normal (ps);
      const LineSurface* ls = dynamic_cast <const LineSurface*> (other);
      if (ls != NULL) return normal (ls);
      const FlatSurface* fs = dynamic_cast <const FlatSurface*> (other);
      if (fs != NULL) return normal (fs);
      assert (false && "Unkown surface primitive");
    }

    AbstractSurface::Point_t PointSurface::normal (const PointSurface* other) const {
      Point_t n = center() - pointToLocalFrame ( 
          other->pointToGlobalFrame (other->center())
          );
      n.normalize();
      return n;
    }
    AbstractSurface::Point_t PointSurface::normal (const LineSurface* other) const {
      return pointToLocalFrame ( 
          other->pointToGlobalFrame (other->normal (this))
          );
    }
    AbstractSurface::Point_t PointSurface::normal (const FlatSurface* other) const {
      return pointToLocalFrame ( 
          other->pointToGlobalFrame (other->normal (this))
          );
    }

    AbstractSurface::Point_t LineSurface::normal (const PointSurface* other) const {
      Point_t cPointInThis = pointToLocalFrame (
          other->pointToGlobalFrame (other->center())
          );
      Point_t n = project (cPointInThis) - cPointInThis;
      n.normalize();
      return n;
    }
    AbstractSurface::Point_t LineSurface::normal (const LineSurface* other) const {
      // The two lines are assumed to be parallel
      Point_t cPointInThis = pointToLocalFrame (
          other->pointToGlobalFrame (other->center())
          );
      Point_t n = project (cPointInThis) - cPointInThis;
      n.normalize();
      return n;
    }
    AbstractSurface::Point_t LineSurface::normal (const FlatSurface* other) const {
      return pointToLocalFrame ( 
          other->pointToGlobalFrame (other->normal (this))
          );
    }

    DifferentiableFunctionPtr_t nullVelocityConstraint (
        const FlatSurface& s1, const FlatSurface& s2, const DevicePtr_t device)
    {
      fcl::Matrix3f R1 (s1.tangentX(), s1.tangentY(), s1.normal ());
      fcl::Transform3f ref1 (R1, s1.center());

      fcl::Matrix3f R2 (s2.tangentX(), s2.tangentY(), s2.normal ());
      fcl::Transform3f ref2 (R2, s2.center());

      if (s1.joint () == NULL) {
        assert (s2.joint () != NULL);
        return Transformation::create (
            "T", device, s2.joint(), ref2,
            boost::assign::list_of (false)(false)(true)(true)(true)(false)
            );
      } else if (s2.joint () == NULL) {
        return Transformation::create (
            "RT", device, s1.joint(), ref1,
            boost::assign::list_of (false)(false)(true)(true)(true)(false)
            );
      } else {
        return RelativeTransformation::create (
            "RT", device, s1.joint(), s2.joint(), ref1, ref2,
            boost::assign::list_of (false)(false)(true)(true)(true)(false)
            );
      }
    }

    DifferentiableFunctionPtr_t nullVelocityConstraint (
        const PointSurface& s1, const FlatSurface& s2, const DevicePtr_t device)
    {
      fcl::Transform3f ref1 (s1.center());

      fcl::Matrix3f R2 (s2.tangentX(), s2.tangentY(), s2.normal ());
      fcl::Transform3f ref2 (R2, s2.center());

      if (s1.joint () == NULL) {
        assert (s2.joint () != NULL);
        return Transformation::create (
            "T", device, s2.joint(), ref2,
            boost::assign::list_of (false)(false)(true)(false)(false)(false)
            );
      } else if (s2.joint () == NULL) {
        return Transformation::create (
            "RT", device, s1.joint(), ref1,
            boost::assign::list_of (false)(false)(true)(false)(false)(false)
            );
      } else {
        return RelativeTransformation::create (
            "RT", device, s1.joint(), s2.joint(), ref1, ref2,
            boost::assign::list_of (false)(false)(true)(false)(false)(false)
            );
      }
    }

    DifferentiableFunctionPtr_t noSlideConstraint (
        const FlatSurface& s1, const FlatSurface& s2, const DevicePtr_t device)
    {
      fcl::Matrix3f R1 (s1.tangentX(), s1.tangentY(), s1.normal ());
      fcl::Transform3f ref1 (R1, s1.center());

      fcl::Matrix3f R2 (s2.tangentX(), s2.tangentY(), s2.normal ());
      fcl::Transform3f ref2 (R2, s2.center());

      if (s1.joint () == NULL) {
        assert (s2.joint () != NULL);
        return Transformation::create (
            "T", device, s2.joint(), ref2,
            boost::assign::list_of (true)(true)(false)(false)(false)(true)
            );
      } else if (s2.joint () == NULL) {
        return Transformation::create (
            "RT", device, s1.joint(), ref1,
            boost::assign::list_of (true)(true)(false)(false)(false)(true)
            );
      } else {
        return RelativeTransformation::create (
            "RT", device, s1.joint(), s2.joint(), ref1, ref2,
            boost::assign::list_of (true)(true)(false)(false)(false)(true)
            );
      }
    }

    DifferentiableFunctionPtr_t noSlideConstraint (
        const PointSurface& s1, const FlatSurface& s2, const DevicePtr_t device)
    {
      fcl::Transform3f ref1 (s1.center());

      fcl::Matrix3f R2 (s2.tangentX(), s2.tangentY(), s2.normal ());
      fcl::Transform3f ref2 (R2, s2.center());

      if (s1.joint () == NULL) {
        assert (s2.joint () != NULL);
        return Transformation::create (
            "T", device, s2.joint(), ref2,
            boost::assign::list_of (true)(true)(false)(false)(false)(false)
            );
      } else if (s2.joint () == NULL) {
        return Transformation::create (
            "RT", device, s1.joint(), ref1,
            boost::assign::list_of (true)(true)(false)(false)(false)(false)
            );
      } else {
        return RelativeTransformation::create (
            "RT", device, s1.joint(), s2.joint(), ref1, ref2,
            boost::assign::list_of (true)(true)(false)(false)(false)(false)
            );
      }
    }

    DifferentiableFunctionPtr_t contactConstraint (
        const FlatSurface& s1, const FlatSurface& s2, const DevicePtr_t device)
    {
    }

    DifferentiableFunctionPtr_t contactConstraint (
        const PointSurface& s1, const FlatSurface& s2, const DevicePtr_t device)
    {
    }
  } // namespace constraints
} // namespace hpp
