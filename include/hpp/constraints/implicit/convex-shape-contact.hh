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

#ifndef HPP_CONSTRAINTS_IMPLICIT_CONVEX_SHAPE_CONTACT_HH
#define HPP_CONSTRAINTS_IMPLICIT_CONVEX_SHAPE_CONTACT_HH

#include <hpp/constraints/implicit.hh>

namespace hpp {
  namespace constraints {
    namespace implicit {
      /// \addtogroup constraints
      /// \{

      /// Contact constraint between sets of convex polygons
      ///
      /// This constraint is the implicit expression of the parameterizable
      /// constraint of contact between a polygon among a set
      /// \f$(f_j)_{j\in J}\f$, and a polygon among another set
      /// \f$(o_i)_{i\in I}\f$.
      ///
      /// The left hand side of the constraint is of type
      /// ConvexShapeContactHold.
      class HPP_CONSTRAINTS_DLLAPI ConvexShapeContact : public virtual Implicit
      {
      public:
        /// Copy object and return shared pointer to copy
        virtual ImplicitPtr_t copy () const;
        /// Create instance and return shared pointer
        /// \param name of the constraint,
        /// \param robot the robot holding the floor and object contact
        ///        surfaces,
        /// \param floorSurfaces set of surfaces on which object are placed,
        /// \param objectSurfaces set of surfaces that come in contact with
        ///        floor surfaces.
        static ConvexShapeContactPtr_t create
          (const std::string& name, DevicePtr_t robot,
           const JointAndShapes_t& floorSurfaces,
           const JointAndShapes_t& objectSurfaces);

        static ConvexShapeContactPtr_t createCopy
          (const ConvexShapeContactPtr_t& other);

      protected:
        /// Constructor
        /// \param name of the constraint,
        /// \param robot the robot holding the floor and object contact
        ///        surfaces,
        /// \param floorSurfaces set of surfaces on which object are placed,
        /// \param objectSurfaces set of surfaces that come in contact with
        ///        floor surfaces.
        ConvexShapeContact
          (const std::string& name, DevicePtr_t robot,
           const JointAndShapes_t& floorSurfaces,
           const JointAndShapes_t& objectSurfaces);
        /// Copy constructor
        ConvexShapeContact(const ConvexShapeContact& other);
        /// Store shared pointer to itself
        void init (ConvexShapeContactWkPtr_t weak);
      private:
        ConvexShapeContactWkPtr_t weak_;
      }; // class ConvexShapeContact
      /// \}
    } // namespace implicit
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_IMPLICIT_CONVEX_SHAPE_CONTACT_HH
