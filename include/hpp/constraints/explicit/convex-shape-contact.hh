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

#ifndef HPP_CONSTRAINTS_EXPLICIT_CONVEX_SHAPE_CONTACT_HH
#define HPP_CONSTRAINTS_EXPLICIT_CONVEX_SHAPE_CONTACT_HH

#include <hpp/constraints/implicit/convex-shape-contact.hh>
#include <hpp/constraints/explicit.hh>

namespace hpp {
  namespace constraints {
    namespace explicit_ {
      typedef implicit::ConvexShapeContact Parent;
      class HPP_CONSTRAINTS_DLLAPI ConvexShapeContact :
        public Explicit, public Parent
      {
      public:
        typedef std::tuple<ImplicitPtr_t, ImplicitPtr_t, ExplicitPtr_t>
          Constraints_t;
        /// Copy object and return shared pointer to copy
        virtual ImplicitPtr_t copy () const;
        /// Create instance and return shared pointer
        /// \param robot the robot holding the floor and object contact
        ///        surfaces,
        /// \param floorSurfaces set of surfaces on which object are placed,
        /// \param objectSurfaces set of surfaces that come in contact with
        ///        floor surfaces.
        /// \param margin distance between surfaces at contact to avoid
        ///        collision.
        static ConvexShapeContactPtr_t create
          (const std::string& name, DevicePtr_t robot,
           const JointAndShapes_t& floorSurfaces,
           const JointAndShapes_t& objectSurfaces,
           const value_type& margin);
        /// Create 3 constraints
        /// \param robot the robot holding the floor and object contact
        ///        surfaces,
        /// \param floorSurfaces set of surfaces on which object are placed,
        /// \param objectSurfaces set of surfaces that come in contact with
        ///        floor surfaces,
        /// \param margin distance between surfaces at contact to avoid
        ///        collision.
        /// \return a tuple containing
        ///         \li a contact constraint (see function ConvexShapeContact),
        ///         \li the complement constraint (see function
        ///             ConvexShapeContactComplement),
        ///         \li an explicit constraint corresponding to the combination
        ///             of the former constraints.
        static Constraints_t createConstraintAndComplement
          (const std::string& name, DevicePtr_t robot,
           const JointAndShapes_t& floorSurfaces,
           const JointAndShapes_t& objectSurfaces,
           const value_type& margin);

        /// Create copy and return shared pointer
        static ConvexShapeContactPtr_t createCopy
          (const ConvexShapeContactPtr_t& other);
        /// Compute the object pose with respect
        /// \param qin input configuration variables,
        /// \param rhs right hand side of constraint
        ///
        /// The right hand side contains the following information:
        /// \li which floor surface is in contact with with object surface,
        /// \li the relative position of those surfaces.
        virtual void outputValue(LiegroupElementRef result, vectorIn_t qin,
                                 vectorIn_t rhs) const;

        /// Compute Jacobian of output value
        ///
        /// \f{eqnarray*}
        /// J &=& \frac{\partial}{\partial\mathbf{q}_{in}}
        ///       \left(f(\mathbf{q}_{in}) + rhs\right).
        /// \f}
        /// \param qin vector of input variables,
        /// \param f_value \f$f(\mathbf{q}_{in})\f$ provided to avoid
        ///        recomputation,
        /// \param rhs right hand side (of implicit formulation).
        virtual void jacobianOutputValue
          (vectorIn_t qin, LiegroupElementConstRef f_value,
           vectorIn_t rhs, matrixOut_t jacobian) const;
      protected:
        /// Constructor
        /// \param robot the robot holding the floor and object contact
        ///        surfaces,
        /// \param floorSurfaces set of surfaces on which object are placed,
        /// \param objectSurfaces set of surfaces that come in contact with
        ///        floor surfaces.
        /// \param margin distance between surfaces at contact to avoid
        ///        collision.
        ConvexShapeContact
          (const std::string& name, DevicePtr_t robot,
           const JointAndShapes_t& floorSurfaces,
           const JointAndShapes_t& objectSurfaces,
           const value_type& margin);
        /// Store weak pointer to itself
        void init (ConvexShapeContactWkPtr_t weak);
      private:
        // Store a vector of explicit relative transforms corresponding to
        // each pair (floor surface, object surface)
        std::vector<RelativePosePtr_t> pose_;
        // Number of floor surfaces
        std::size_t nFloor_;
        // shared pointer to itself
        ConvexShapeContactWkPtr_t weak_;
      }; // class ConvexShapeContact
    } // namespace explicit_
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_EXPLICIT_CONVEX_SHAPE_CONTACT_HH
