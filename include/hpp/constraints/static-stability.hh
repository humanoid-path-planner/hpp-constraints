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

#ifndef HPP_CONSTRAINTS_STATIC_STABILITY_HH
# define HPP_CONSTRAINTS_STATIC_STABILITY_HH

# include <vector>
# include <hpp/constraints/differentiable-function.hh>
# include <hpp/fcl/math/transform.h>
# include <hpp/fcl/shape/geometric_shapes.h>
# include <hpp/constraints/tools.hh>
# include <hpp/constraints/convex-shape.hh>

# include "hpp/constraints/fwd.hh"
# include "hpp/constraints/config.hh"
# include "hpp/constraints/deprecated.hh"

namespace hpp {
  namespace constraints {
    class HPP_CONSTRAINTS_DLLAPI Triangle {
      public:
        /// Represent a triangle.
        Triangle (const fcl::Vec3f& p0, const fcl::Vec3f& p1, const fcl::Vec3f& p2,
            JointPtr_t joint = NULL):
          P0_ (p0), P1_ (p1), P2_ (p2), joint_ (joint),
          p0_ (p0), p1_ (p1), p2_ (p2)
        {
          init ();
          recompute ();
        }

        Triangle (const fcl::TriangleP& t, const JointPtr_t& joint = NULL):
          P0_ (t.a), P1_ (t.b), P2_ (t.c), joint_ (joint),
          p0_ (t.a), p1_ (t.b), p2_ (t.c)
        {
          init ();
          recompute ();
        }

        // Copy constructor
        Triangle (const Triangle& t) :
          P0_ (t.P0_), P1_ (t.P1_), P2_ (t.P2_), joint_ (t.joint_),
          p0_ (t.P0_), p1_ (t.P1_), p2_ (t.P2_)
        {
          init ();
          recompute ();
        }

        void updateToCurrentTransform () const
        {
          if (joint_ == NULL) return;
          const Transform3f& M = joint_->currentTransformation ();
          p0_ = M.transform (P0_);
          p1_ = M.transform (P1_);
          p2_ = M.transform (P2_);
          recompute ();
        }

        /// Intersection with a line defined by a point and a vector.
        inline fcl::Vec3f intersection (const fcl::Vec3f& A, const fcl::Vec3f& u) const {
          assert (std::abs (n_.dot (u)) > 1e-8);
          return A + u * (n_.dot(c_ - A)) / n_.dot (u);
        }

        /// Check whether the intersection of the line defined by A and u
        /// onto the plane containing the triangle is inside the triangle.
        inline bool isInside (const fcl::Vec3f& A, const fcl::Vec3f& u) const {
          fcl::Vec3f Ap = intersection (A, u);
          return isInside (Ap);
        }
        inline bool isInside (const fcl::Vec3f& Ap) const {
          if (n2_.dot (Ap-p0_) > 0) return false;
          if (n0_.dot (Ap-p1_) > 0) return false;
          if (n1_.dot (Ap-p2_) > 0) return false;
          return true;
        }

        /// Return the shortest distance from a point to the triangle.
        /// A negative value means the point is inside the triangle.
        /// \param A a point already in the plane containing the triangle.
        inline value_type distance (const fcl::Vec3f& A) const {
          value_type d0 = dist (A-p0_, p1_-p0_, n2_),
                     d1 = dist (A-p1_, p2_-p1_, n0_),
                     d2 = dist (A-p2_, p0_-p2_, n1_);
          // We should return the shortest positive distance
          bool d0p = d0 >= 0;
          bool d1p = d1 >= 0;
          bool d2p = d2 >= 0;
          if (d0p) {
            if (d2 > d1) return (d2 > d0)?d2:d0;
            else return (d1 > d0)?d1:d0;
          } else if (d1p) {
            if (d2 > d0) return (d2 > d1)?d2:d1;
            else return (d0 > d1)?d0:d1;
          } else if (d2p) {
            if (d1 > d0) return (d1 > d2)?d1:d2;
            else return (d0 > d2)?d0:d2;
          } else {
            if (d0 < d1) return (d2 > d1)?d2:d1;
            else return (d0 > d2)?d0:d2;
          }
        }

        inline const fcl::Vec3f& planeXaxis () const { return   n0_; }
        inline const fcl::Vec3f& planeYaxis () const { return nxn0_; }
        inline const fcl::Vec3f& normal () const { return n_; }
        inline const fcl::Vec3f& center () const { return c_; }
      /// Return inverse of position of triangle local frame
      ///
      /// Triangle local frame is defined by
      /// \li center is barycenter of 3 vertices,
      /// \li x-axis is normal to triangle plane,
      /// \li y-axis is aligned with \f$(p2 - p1)\times \vec{x}\f$
      /// \li z-axis is equal to \f$\vec{x}\times\vec{y}\f$.
        inline const fcl::Transform3f& inversePosition () const { return M_; }

      /// Write triangle in a stream
      std::ostream& print (std::ostream& os) const
      {
	os << "triangle: (" << p0_ << "," << p1_ << "," << p2_ << ")";
	return os;
      }

        /// The position in the joint frame and the joint
        fcl::Vec3f P0_, P1_, P2_, C_;
        JointPtr_t joint_;

      private:
        /// Return the distance between the point A and the segment
        /// [P, v] oriented by u.
        /// w = PA.
        inline value_type dist (const fcl::Vec3f& w, const fcl::Vec3f& v, const fcl::Vec3f& u) const {
          value_type c1, c2;
          c1 = v.dot (w);
          if (c1 <= 0)
            return (u.dot (w) > 0)?(w.norm()):(- w.norm());
          c2 = v.norm ();
          if (c2 <= c1)
            return (u.dot (w) > 0)?((w-v).norm()):(-(w-v).norm());
          return u.dot (w);
        }

        void init ()
        {
          C_ = ( P0_ + P1_ + P2_ ) / 3;
        }

        void recompute () const
        {
          n_ = (p1_ - p0_).cross (p2_ - p0_);
          assert (!n_.isZero ());
          n_.normalize ();
          c_ = ( p0_ + p1_ + p2_ ) / 3;
          n0_ = (p2_ - p1_).cross (n_); n0_.normalize ();
          n1_ = (p0_ - p2_).cross (n_); n1_.normalize ();
          n2_ = (p1_ - p0_).cross (n_); n2_.normalize ();
          nxn0_ = n_.cross (n0_);
          M_ = fcl::Transform3f (fcl::Matrix3f (n_, n0_, nxn0_));
          M_.setTranslation (- (M_.getRotation () * c_));
          for (size_t i = 0; i < 3; i++) assert (M_.getRotation () (0, i) == n_[i]);
          for (size_t i = 0; i < 3; i++) assert (M_.getRotation () (1, i) == n0_[i]);
          for (size_t i = 0; i < 3; i++) assert (M_.getRotation () (2, i) == nxn0_[i]);
        }

        /// The positions and vectors in the global frame
        mutable fcl::Vec3f p0_, p1_, p2_, n_, c_;
        /// n_i is the vector of norm 1 perpendicular to
        /// P_{i+1}P_i and n_.
        mutable fcl::Vec3f n0_, n1_, n2_, nxn0_;
        mutable fcl::Transform3f M_;
    };
    std::ostream& operator<< (std::ostream& os, const Triangle& t);

    /// \addtogroup constraints
    /// \{

    class HPP_CONSTRAINTS_DLLAPI StaticStabilityGravity :
      public DifferentiableFunction {
      public:
      friend class StaticStabilityGravityComplement;
        /// Constructor
        /// \param robot the robot the constraints is applied to,
        /// \param com COM of the object in the joint frame.
        StaticStabilityGravity (const std::string& name,
				const DevicePtr_t& robot);

        static StaticStabilityGravityPtr_t create (
            const std::string& name,
            const DevicePtr_t& robot);

        static StaticStabilityGravityPtr_t create (
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

        void addObject (const ConvexShape& t);

        void addFloor (const ConvexShape& t);

      private:
        enum ContactType {
          POINT_ON_PLANE,
          LINE_ON_PLANE,
          PLANE_ON_PLANE
        };

        void impl_compute (vectorOut_t result, ConfigurationIn_t argument) const;

        void impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const;
        void computeInternalJacobian (ConfigurationIn_t argument) const;

        void selectConvexShapes () const;
        ContactType contactType (const ConvexShape& object,
            const ConvexShape& floor) const;

        DevicePtr_t robot_;
        RelativeTransformationPtr_t relativeTransformation_;

        typedef std::vector <ConvexShape> ConvexShapes_t;
        ConvexShapes_t objectConvexShapes_;
        ConvexShapes_t floorConvexShapes_;

        mutable bool isInside_;
        mutable ContactType contactType_;
        mutable vector_t result_;
        mutable matrix_t jacobian_;
    };

    /// Complement to full transformation constraint of StaticStabilityGravity
    class HPP_CONSTRAINTS_DLLAPI StaticStabilityGravityComplement :
      public DifferentiableFunction
    {
    public:
      /// Create a pair of constraints
      ///
      /// The pair contains two complementary constraints to be used for
      /// manipulation applications.
      /// \param name name of the static stability constraint,
      /// \param constraintName name of the complement constraint,
      /// \param name of the robot.
      static std::pair <StaticStabilityGravityPtr_t,
			StaticStabilityGravityComplementPtr_t >
	createPair (const std::string& name, const std::string& complementName,
		    const DevicePtr_t& robot);

    protected:
      /// Constructor
      /// \param name name of the static stability constraint,
      /// \param constraintName name of the complement constraint,
      /// \param name of the robot.
      StaticStabilityGravityComplement (const std::string& name,
					const std::string& complementName,
					const DevicePtr_t& robot);


    private:
      void impl_compute (vectorOut_t result, ConfigurationIn_t argument) const;

      void impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument)
	const;

      StaticStabilityGravityPtr_t sibling_;
    }; // class StaticStabilityGravityComplement

    /// \}

    /// \addtogroup constraints
    /// \{

    class HPP_CONSTRAINTS_DLLAPI StaticStability :
      public DifferentiableFunction {
      public:
        static const value_type G;
        static const Eigen::Matrix <value_type, 6, 1> Gravity;

        struct Contact_t {
          JointPtr_t joint1, joint2;
          vector3_t point1, point2;
          vector3_t normal1, normal2;
        };
        typedef std::vector <Contact_t> Contacts_t;

        /// Constructor
        /// \param robot the robot the constraints is applied to,
        /// \param joint the joint to which the object is attached,
        /// \param com COM of the object in the joint frame.
        StaticStability (const std::string& name, const DevicePtr_t& robot,
            const Contacts_t& contacts,
            const CenterOfMassComputationPtr_t& com);

        static StaticStabilityPtr_t create (
            const std::string& name,
            const DevicePtr_t& robot,
            const Contacts_t& contacts,
            const CenterOfMassComputationPtr_t& com);

        static StaticStabilityPtr_t create (
            const DevicePtr_t& robot,
            const Contacts_t& contacts,
            const CenterOfMassComputationPtr_t& com);

        MatrixOfExpressions<>& phi () {
          return phi_;
        }

      private:
        void impl_compute (vectorOut_t result, ConfigurationIn_t argument) const;

        void impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const;

        static void findBoundIndex (vectorIn_t u, vectorIn_t v, 
            value_type& lambdaMin, size_type* iMin,
            value_type& lambdaMax, size_type* iMax);

        /// Return false if uMinus.isZero(), i which case v also zero (not computed).
        bool computeUminusAndV (vectorIn_t u, vectorOut_t uMinus,
            vectorOut_t v) const;

        void computeVDot (vectorIn_t uMinus, vectorIn_t S,
            matrixIn_t uDot, matrixOut_t uMinusDot, matrixOut_t vDot) const;

        void computeLambdaDot (vectorIn_t u, vectorIn_t v, const std::size_t i0,
            matrixIn_t uDot, matrixIn_t vDot, vectorOut_t lambdaDot) const;

        DevicePtr_t robot_;
        Contacts_t contacts_;
        CenterOfMassComputationPtr_t com_;

        typedef MatrixOfExpressions<eigen::vector3_t, JacobianMatrix> MoE_t;

        mutable MoE_t phi_;
        mutable vector_t u_, uMinus_, v_;
        mutable matrix_t uDot_, uMinusDot_, vDot_;
        mutable vector_t lambdaDot_; 
    };
    /// \}
  } // namespace constraints
} // namespace hpp

#endif //  HPP_CONSTRAINTS_STATIC_STABILITY_HH
