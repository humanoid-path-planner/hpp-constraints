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

#include <boost/assign/list_of.hpp>
#include <hpp/constraints/convex-shape-contact.hh>
#include <hpp/constraints/implicit/convex-shape-contact.hh>

namespace hpp {
  namespace constraints {
    namespace implicit {
      using boost::assign::list_of;
      ImplicitPtr_t ConvexShapeContact::copy () const
      {
        return createCopy (weak_.lock ());
      }

      ConvexShapeContactPtr_t ConvexShapeContact::create
      (const std::string& name, DevicePtr_t robot,
       const JointAndShapes_t& floorSurfaces,
       const JointAndShapes_t& objectSurfaces)
      {
        ConvexShapeContact* ptr (new ConvexShapeContact
                                 (name, robot, floorSurfaces, objectSurfaces));
        ConvexShapeContactPtr_t shPtr (ptr);
        ConvexShapeContactWkPtr_t wkPtr (shPtr);
        ptr->init (shPtr);
        return shPtr;
      }

      ConvexShapeContactPtr_t ConvexShapeContact::createCopy
      (const ConvexShapeContactPtr_t& other)
      {
        ConvexShapeContact* ptr (new ConvexShapeContact (*other));
        ConvexShapeContactPtr_t shPtr (ptr);
        ConvexShapeContactWkPtr_t wkPtr (shPtr);
        ptr->init (shPtr);
        return shPtr;
      }

      ConvexShapeContact::ConvexShapeContact
      (const std::string& name, DevicePtr_t robot,
       const JointAndShapes_t& floorSurfaces,
       const JointAndShapes_t& objectSurfaces) :
        Implicit (ConvexShapeContactHold::create(name, robot, floorSurfaces,
                                                 objectSurfaces),
                  list_of(EqualToZero)(EqualToZero)(EqualToZero)(EqualToZero)
                  (EqualToZero)(Equality)(Equality)(Equality))
      {
      }

      ConvexShapeContact::ConvexShapeContact (const ConvexShapeContact& other) :
        Implicit (other)
      {
      }

      void ConvexShapeContact::init (ConvexShapeContactWkPtr_t weak)
      {
        Implicit::init (weak);
        weak_ = weak;
      }
    } // namespace implicit
  } // namespace constraints
} // namespace hpp
