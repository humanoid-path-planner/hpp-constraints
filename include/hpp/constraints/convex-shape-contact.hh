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

        Two set of convex shapes can be given to this class:
        \li a set of object contact surfaces, \f$ (o_i)_{i \in I } \f$, which can be in contact with the environment,
        \li a set of floor contact surfaces, \f$ (f_j)_{j \in J } \f$, which can support objects.

        The distance \f$ d_{i,j} = d (o_i, f_j) \f$ between object surface
        \f$o_i\f$ and environment surface \f$ f_j \f$ is defined by:
        \f{equation*}
           d(i,j)^2 =
             \left\lbrace \begin{array}{cl}
               d_{\parallel}^2 + d_{\perp}^2 &, \text{ if } d_{\parallel} > 0 \\
               d_{\perp}^2                   &, \text{ otherwise}
             \end{array} \right.
        \f}
        where
        \li \f$P (C_{o_i}, f_j)\f$ is the projection of the center \f$o_i\f$ onto the plane containing \f$ f_j \f$,
        \li \f$\textbf{n}_{f_j}\f$ is the normal of \f$ f_j \f$,
        \li \f$d_{\parallel} = d(f_j, P (C_{o_i}, f_j))\f$ is the distance returned by ConvexShape::distance,
        \li \f$d_{\perp} = \textbf{n}_{f_j}.C_{f_j}P(C_{o_i}, f_j)\f$ is the distance along the normal of \f$ f_j \f$,

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

        /// Constructor
        /// \param name name of the ConvexShapeContact constraint,
        /// \param robot the robot the constraints is applied to,
        ConvexShapeContact (const std::string& name,
				const DevicePtr_t& robot);

        static ConvexShapeContactPtr_t create (
            const std::string& name,
            const DevicePtr_t& robot);

        static ConvexShapeContactPtr_t create (
            const DevicePtr_t& robot);

        /// Use addObject(const ConvexShape&) instead.
        /// Add a triangle to the object contact surface
        /// \param t triangle,
        /// \param joint Joint to which the triangle is attached.
        void addObjectTriangle (const fcl::TriangleP& t,
				const JointPtr_t& joint)
          HPP_CONSTRAINTS_DEPRECATED;

        /// Use addFloor(const ConvexShape&) instead.
        /// Add a triangle to the floor contact surface
        /// \param t triangle,
        /// joint Joint to which the triangle is attached if the contact surface
        ///       belongs to a robot.
        void addFloorTriangle (const fcl::TriangleP& t,
			       const JointPtr_t& joint)
          HPP_CONSTRAINTS_DEPRECATED;

        /// Add a ConvexShape as an object.
        void addObject (const ConvexShape& t);

        /// Add a ConvexShape as a floor.
        ///
        /// The convex shape will be reverted using ConvexShape::reverse
        /// so that the normal points inside the floor object.
        void addFloor (const ConvexShape& t);

        /// Set the normal margin, i.e. the desired distance between matching
        /// object and nd floor shapes.
        /// Default to 0
        void setNormalMargin (const value_type& margin);

        /// Compute the contact points
        std::vector <ForceData> computeContactPoints (ConfigurationIn_t q,
            const value_type& normalMargin) const;

      /// Display object in a stream
      std::ostream& print (std::ostream& o) const;

      private:
        void impl_compute (LiegroupElementRef result, ConfigurationIn_t argument)
          const;
        void computeInternalValue (const ConfigurationIn_t& argument,
            bool& isInside, ContactType& type, vector6_t& value) const;

        void impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const;
        void computeInternalJacobian (const ConfigurationIn_t& argument,
            bool& isInside, ContactType& type, matrix_t& jacobian) const;

        typedef std::vector <ConvexShape> ConvexShapes_t;
        /// \return true if the contact is created.
        bool selectConvexShapes (const pinocchio::DeviceData& data,
            ConvexShapes_t::const_iterator& object,
            ConvexShapes_t::const_iterator& floor) const;
        ContactType contactType (const ConvexShape& object,
            const ConvexShape& floor) const;

        DevicePtr_t robot_;
        mutable GenericTransformationModel<true> relativeTransformationModel_;

        ConvexShapes_t objectConvexShapes_;
        ConvexShapes_t floorConvexShapes_;

        value_type normalMargin_;
    };

    /** Complement to full transformation constraint of ConvexShapeContact

        The value returned by this class is:

        | Contact type   | Inside   | Outside |
        | -------------- | -------- | ------- |
        | ConvexShapeContact::POINT_ON_PLANE (Unsupported) | \f$(y,z,rx)\f$ | \f$(0,0,rx)\f$ |
        | ConvexShapeContact::LINE_ON_PLANE (Unsupported)  | \f$(y,z,rx)\f$ | \f$(0,0,rx)\f$  |
        | ConvexShapeContact::PLANE_ON_PLANE | \f$(y,z,rx)\f$ | \f$(0,0,rx)\f$ |

        See ConvexShapeContact
     **/
    class HPP_CONSTRAINTS_DLLAPI ConvexShapeContactComplement :
      public DifferentiableFunction
    {
    public:
      /// Create a pair of constraints
      ///
      /// The pair contains two complementary constraints to be used for
      /// manipulation applications.
      /// \param name name of the ConvexShapeContact constraint,
      /// \param complementName name of the complement constraint,
      /// \param robot
      static std::pair <ConvexShapeContactPtr_t,
			ConvexShapeContactComplementPtr_t >
	createPair (const std::string& name, const std::string& complementName,
		    const DevicePtr_t& robot);

    protected:
      /// Constructor
      /// \param name name of the ConvexShapeContact constraint,
      /// \param complementName name of the complement constraint,
      /// \param robot
      ConvexShapeContactComplement (const std::string& name,
					const std::string& complementName,
					const DevicePtr_t& robot);


    private:
      void impl_compute (LiegroupElementRef result, ConfigurationIn_t argument) const;

      void impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument)
	const;

      ConvexShapeContactPtr_t sibling_;
    }; // class ConvexShapeContactComplement
    /// \}
  } // namespace constraints
} // namespace hpp

#endif //  HPP_CONSTRAINTS_CONVEX_SHAPE_CONTACT_HH
