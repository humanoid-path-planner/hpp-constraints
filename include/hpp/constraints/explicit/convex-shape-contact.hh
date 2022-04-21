// Copyright (c) 2020, Airbus SAS and CNRS
// Authors: Florent Lamiraux
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

#ifndef HPP_CONSTRAINTS_EXPLICIT_CONVEX_SHAPE_CONTACT_HH
#define HPP_CONSTRAINTS_EXPLICIT_CONVEX_SHAPE_CONTACT_HH

#include <hpp/constraints/explicit.hh>
#include <tuple>

namespace hpp {
  namespace constraints {
    namespace explicit_ {
      class HPP_CONSTRAINTS_DLLAPI ConvexShapeContact : public Explicit
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
        /// \param rhs right hand side of constraint (of implicit formulation).
        ///
        /// The right hand side contains the following information:
        /// \li which floor surface is in contact with with object surface,
        /// \li the relative position of those surfaces.
        virtual void outputValue(LiegroupElementRef result, vectorIn_t qin,
                                 LiegroupElementConstRef rhs) const;

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
           LiegroupElementConstRef rhs, matrixOut_t jacobian) const;
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
        // Number of object surfaces
        std::size_t nObjSurf_;
        // shared pointer to itself
        ConvexShapeContactWkPtr_t weak_;
      }; // class ConvexShapeContact
    } // namespace explicit_
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_EXPLICIT_CONVEX_SHAPE_CONTACT_HH
