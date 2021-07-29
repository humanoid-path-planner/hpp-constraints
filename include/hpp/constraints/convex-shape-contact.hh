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

#ifndef HPP_CONSTRAINTS_CONVEX_SHAPE_CONTACT_HH
# define HPP_CONSTRAINTS_CONVEX_SHAPE_CONTACT_HH

# include <vector>

# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/config.hh>
# include <hpp/constraints/deprecated.hh>
# include <hpp/constraints/generic-transformation.hh>
# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/convex-shape.hh>

namespace hpp {
  namespace constraints {

    /// \addtogroup constraints
    /// \{

    /** The function returns a relative transformation between the two "closest"
        convex shapes it contains.

        Twos set of convex shapes can be given to this class:
        \li a set of object contact surfaces, \f$ (o_i)_{i \in I } \f$, which can be in contact with the environment,
        \li a set of floor contact surfaces, \f$ (f_j)_{j \in J } \f$, which can support objects.

        The distance \f$ d (f_j, o_i) \f$ between object surface
        \f$o_i\f$ and environment surface \f$ f_j \f$ is defined by:
        \f{equation*}
           d (f_j, o_i)^2 =
             \left\lbrace \begin{array}{cl}
               d_{\parallel}^2 + d_{\perp}^2 &, \text{ if } d_{\parallel} > 0 \\
               d_{\perp}^2                   &, \text{ otherwise}
             \end{array} \right.
        \f}
        where
        \li \f$P (C_{o_i}, f_j)\f$ is the projection of the center \f$C_{o_i}\f$ of \f$o_i\f$ onto the plane containing \f$ f_j \f$,
        \li \f$\textbf{n}_{f_j}\f$ is the normal of \f$ f_j \f$,
        \li \f$d_{\parallel} = d(f_j, P (C_{o_i}, f_j))\f$ is the distance returned by ConvexShapeData::distance,
        \li \f$d_{\perp} = \textbf{n}_{f_j}.\vec{C_{f_j}C_{o_i}}\f$ is the distance along the normal of \f$ f_j \f$,

        \image html convex-shape-contact.svg

        The function first selects the pair \f$(o_i,f_j)\f$ with shortest distance.
        \f$o_i\f$ is \em inside \f$f_j\f$ if \f$d(i,j) < 0\f$.
        It returns a value that depends on the contact types:


        | Contact type   | Inside   | Outside |
        | -------------- | -------- | ------- |
        | ConvexShapeContact::POINT_ON_PLANE | \f$(x+m,0,0,0,0)\f$ | \f$(x+m,y,z,0,0)\f$ |
        | ConvexShapeContact::LINE_ON_PLANE (Unsupported)  | \f$(x+m,0,0,0,rz)\f$ | \f$(x+m,y,z,0,rz)\f$  |
        | ConvexShapeContact::PLANE_ON_PLANE | \f$(x+m,0,0,ry,rz)\f$ | \f$(x+m,y,z,ry,rz)\f$ |

        where
        \li \f$m\f$ is the normal margin (used to avoid collisions),
        \li \f$x,y,z,rx,ry,rz\f$ represents the output of the RelativeTransformation
            between the element of the pair.

        \sa ConvexShapeContactComplement
    **/
    class HPP_CONSTRAINTS_DLLAPI ConvexShapeContact :
      public DifferentiableFunction {
      public:
        friend class ConvexShapeContactComplement;
        friend class ConvexShapeContactHold;

        /// \cond
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// \endcond

        /// The type of contact between each pair (object shape, floor shape).
        enum ContactType {
          /// The object shape is a single point,
          POINT_ON_PLANE,
          /// The object shape degenerates to a line,
          LINE_ON_PLANE,
          /// The object shape is included in a plane and none of the above case apply.
          PLANE_ON_PLANE
        };

        /// Represents a contact
        /// When supportJoint is NULL, the contact is with the environment.
        /// Otherwise, the contact is between two joints.
        struct ForceData {
          JointPtr_t joint;
          JointPtr_t supportJoint;
          std::vector<vector3_t> points;
          vector3_t normal;
        };

        /// Create instance and return shared pointer
        /// \param name name of the constraint,
        /// \param robot robot that holds the contact surface
        /// \param floorSurfaces, objectSurfaces set of plane polygonal contact
        ///        surfaces.
        static ConvexShapeContactPtr_t create
          (const std::string& name, DevicePtr_t robot,
           const JointAndShapes_t& floorSurfaces,
           const JointAndShapes_t& objectSurfaces);

        static ConvexShapeContactPtr_t create (
            const DevicePtr_t& robot);

        /// Get vector of floor contact surfaces
        const ConvexShapes_t& floorContactSurfaces () const
        {
          return floorConvexShapes_;
        }
        /// Get vector of object contact surfaces
        const ConvexShapes_t& objectContactSurfaces () const
        {
          return objectConvexShapes_;
        }
        /// Get radius \f$M\f$
        /// See class documentation for a definition.
        value_type radius() const
        {
          return M_;
        }
        /// Set the normal margin, i.e. the desired distance between matching
        /// object and nd floor shapes.
        /// Default to 0
        void setNormalMargin (const value_type& margin);

        /// Compute the contact points
        std::vector <ForceData> computeContactPoints (ConfigurationIn_t q,
            const value_type& normalMargin) const;

      /// Display object in a stream
      std::ostream& print (std::ostream& o) const;
    protected:
      /// Constructor
      /// \param name name of the constraint,
      /// \param robot robot that holds the contact surface
      /// \param floorSurfaces, objectSurfaces set of plane polygonal contact
      ///        surfaces.
      ConvexShapeContact (const std::string& name, DevicePtr_t robot,
                          const JointAndShapes_t& floorSurfaces,
                          const JointAndShapes_t& objectSurfaces);

      bool isEqual(const DifferentiableFunction& other) const {
        const ConvexShapeContact& castother = dynamic_cast<const ConvexShapeContact&>(other);
        if (!DifferentiableFunction::isEqual(other))
          return false;
        
        if (robot_ != castother.robot_)
          return false;
        if (objectConvexShapes_.size() != castother.objectConvexShapes_.size())
          return false;
        for (std::size_t i = 0; i < objectConvexShapes_.size(); i++)
          if (floorConvexShapes_[i] != castother.floorConvexShapes_[i])
            return false;
        
        return true;
      }

      private:
        /// Add a ConvexShape as an object.
        void addObject (const ConvexShape& t);

        /// Add a ConvexShape as a floor.
        ///
        /// The convex shape will be reverted using ConvexShape::reverse
        /// so that the normal points inside the floor object.
        void addFloor (const ConvexShape& t);
        void computeRadius();

        void impl_compute (LiegroupElementRef result, ConfigurationIn_t argument)
          const;
        void computeInternalValue (const ConfigurationIn_t& argument,
          bool& isInside, ContactType& type, vector6_t& value,
          std::size_t& iobject, std::size_t& ifloor) const;

        void impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const;
        void computeInternalJacobian (const ConfigurationIn_t& argument,
            bool& isInside, ContactType& type, matrix_t& jacobian) const;

        /// Find floor and object surfaces that are the closest.
        /// \retval iobject, ifloor indices in internal vectors
        ///         objectConvexShapes_ and floorConvexShapes_
        /// \return true if the contact is created.
        bool selectConvexShapes (const pinocchio::DeviceData& data,
                                 std::size_t& iobject, std::size_t& ifloor)
          const;
        ContactType contactType (const ConvexShape& object,
            const ConvexShape& floor) const;

        DevicePtr_t robot_;
        mutable GenericTransformationModel<true> relativeTransformationModel_;

        ConvexShapes_t objectConvexShapes_;
        ConvexShapes_t floorConvexShapes_;

        value_type normalMargin_;
        // upper bound of distance between center of polygon and vectices for
        // all floor polygons.
        value_type M_;
    };

    /** Complement to full transformation constraint of ConvexShapeContact

        The value returned by this class is:

        | Contact type   | Inside   | Outside |
        | -------------- | -------- | ------- |
        | ConvexShapeContact::POINT_ON_PLANE (Unsupported) | \f$(y,z,rx)\f$ | \f$(0,0,rx)\f$ |
        | ConvexShapeContact::LINE_ON_PLANE (Unsupported)  | \f$(y,z,rx)\f$ | \f$(0,0,rx)\f$  |
        | ConvexShapeContact::PLANE_ON_PLANE | \f$(y+2jM,z+2iM,rx)\f$ | \f$(2jM,2iM,rx)\f$ |

        where
        \li \f$M\f$ is an upper bound on the radius of all floor polygons,
        \li \f$i\f$ and \f$j\f$ are the indices of object and floor contact
        surfaces that minimize
        \f{equation}
        d(o_i,f_j)
        \f}

        \sa ConvexShapeContact
     **/
    class HPP_CONSTRAINTS_DLLAPI ConvexShapeContactComplement :
      public DifferentiableFunction
    {
    public:
      friend class ConvexShapeContactHold;
      /// Create a pair of constraints
      /// \param name name of the sibling ConvexShapeContact constraint,
      ///        "/complement" is added to the name of this constraint
      /// \param robot robot that holds the contact surface
      /// \param floorSurfaces, objectSurfaces set of plane polygonal contact
      ///        surfaces.
      static std::pair <ConvexShapeContactPtr_t,
			ConvexShapeContactComplementPtr_t >
	createPair(const std::string& name, DevicePtr_t robot,
                   const JointAndShapes_t& floorSurfaces,
                   const JointAndShapes_t& objectSurfaces);

      /// Compute parameters and right hand side of relative pose
      /// \param rhs right hand side of this constraint,
      /// \retval ifloor, iobject indices of floor and object contact surface,
      /// \retval relativePoseRhs right hand side (implicit representation) of
      ///         the relative pose constraint corresponding to contact of
      ///         surfaces ifloor with iobject.
      void computeRelativePoseRightHandSide
        (LiegroupElementConstRef rhs, std::size_t& ifloor, std::size_t& iobject,
         LiegroupElementRef relativePoseRhs) const;
    protected:
      /// Constructor
      /// \param name name of the sibling ConvexShapeContact constraint,
      ///        "/complement" is added to the name of this constraint
      /// \param robot robot that holds the contact surface
      /// \param floorSurfaces, objectSurfaces set of plane polygonal contact
      ///        surfaces.
      ConvexShapeContactComplement (const std::string& name, DevicePtr_t robot,
                                    const JointAndShapes_t& floorSurfaces,
                                    const JointAndShapes_t& objectSurfaces);

      bool isEqual(const DifferentiableFunction& other) const {
        const ConvexShapeContactComplement& castother = dynamic_cast<const ConvexShapeContactComplement&>(other);
        if (!DifferentiableFunction::isEqual(other))
          return false;
        return (sibling_ != castother.sibling_);
      }
    private:
      void impl_compute (LiegroupElementRef result, ConfigurationIn_t argument) const;

      void impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument)
	const;

      ConvexShapeContactPtr_t sibling_;
    }; // class ConvexShapeContactComplement

    /// Combination of ConvexShapeContact and complement constaints
    ///
    /// Create one instance of
    /// \li ConvexShapeContact,
    /// \li ConvexShapeContactComplement,
    ///
    /// and concatenate their values.
    class HPP_CONSTRAINTS_DLLAPI ConvexShapeContactHold :
      public DifferentiableFunction
    {
    public:
      /// Create instance and return shared pointer
      /// \param constraintName name of the ConvexShapeContact instance,
      /// \param complementName name of the ConvexShapeContactComplement
      ///        instance,
      /// \param holdName name of this DifferentiableFunction instance
      /// \param robot the input space of the function is the robot
      ///        configuration space.
      static ConvexShapeContactHoldPtr_t create
        (const std::string& name, DevicePtr_t robot,
         const JointAndShapes_t& floorSurfaces,
         const JointAndShapes_t& objectSurfaces);

      ConvexShapeContactPtr_t contactConstraint() const
      {
        return constraint_;
      }
      ConvexShapeContactComplementPtr_t complement() const
      {
        return complement_;
      }

    protected:
      /// Constructor
      /// \param constraintName name of the ConvexShapeContact instance,
      /// \param complementName name of the ConvexShapeContactComplement
      ///        instance,
      /// \param holdName name of this DifferentiableFunction instance
      /// \param robot the input space of the function is the robot
      ///        configuration space.
      ConvexShapeContactHold
        (const std::string& name, DevicePtr_t robot,
         const JointAndShapes_t& floorSurfaces,
         const JointAndShapes_t& objectSurfaces);

      virtual void impl_compute(LiegroupElementRef result, vectorIn_t argument)
        const;
      virtual void impl_jacobian(matrixOut_t jacobian, vectorIn_t arg)
        const;

      bool isEqual(const DifferentiableFunction& other) const {
        const ConvexShapeContactHold& castother = dynamic_cast<const ConvexShapeContactHold&>(other);
        if (!DifferentiableFunction::isEqual(other))
          return false;
        
        if (constraint_ != castother.constraint_)
          return false;
        if (complement_ != castother.complement_)
          return false;
        
        return true;
      }
    private:
      ConvexShapeContactPtr_t constraint_;
      ConvexShapeContactComplementPtr_t complement_;
    }; // class ConvexShapeContactHold
    /// \}
  } // namespace constraints
} // namespace hpp

#endif //  HPP_CONSTRAINTS_CONVEX_SHAPE_CONTACT_HH
