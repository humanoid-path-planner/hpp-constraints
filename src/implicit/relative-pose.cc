// Copyright (c) 2018, LAAS-CNRS
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

#include <hpp/constraints/implicit/relative-pose.hh>
#include <hpp/constraints/generic-transformation.hh>

namespace hpp {
  namespace constraints {
    namespace implicit {

      ImplicitPtr_t RelativePose::copy () const
      {
        return createCopy (weak_.lock ());
      }

      RelativePosePtr_t RelativePose::create
      (const std::string& name, const DevicePtr_t& robot,
       const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
       const Transform3f& frame1, const Transform3f& frame2,
       std::vector <bool> mask, ComparisonTypes_t comp, vectorIn_t rhs)
      {
        RelativePose* ptr (new RelativePose
                           (name, robot, joint1, joint2, frame1,
                            frame2, mask, comp, rhs));
        RelativePosePtr_t shPtr (ptr);
        RelativePoseWkPtr_t wkPtr (shPtr);
        ptr->init (shPtr);
        return shPtr;
      }

      RelativePosePtr_t RelativePose::createCopy
      (const RelativePosePtr_t& other)
      {
        RelativePose* ptr (new RelativePose (*other));
        RelativePosePtr_t shPtr (ptr);
        RelativePoseWkPtr_t wkPtr (shPtr);
        ptr->init (shPtr);
        return shPtr;
      }

      RelativePose::RelativePose
      (const std::string& name, const DevicePtr_t& robot,
       const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
       const Transform3f& frame1, const Transform3f& frame2,
       std::vector <bool> mask, ComparisonTypes_t comp, vectorIn_t rhs) :
        Implicit (GenericTransformation
                  <RelativeBit | PositionBit | OrientationBit>::create
                  (name, robot, joint1, joint2, frame1, frame2, mask),
                  comp, rhs), joint1_ (joint1), joint2_ (joint2)
      {
      }

      RelativePose::RelativePose (const RelativePose& other) :
        Implicit (other), joint1_ (other.joint1_), joint2_ (other.joint2_)
      {
      }

      void RelativePose::init (RelativePoseWkPtr_t weak)
      {
        Implicit::init (weak);
        weak_ = weak;
      }
    } // namespace implicit
  } // namespace constraints
} // namespace hpp
