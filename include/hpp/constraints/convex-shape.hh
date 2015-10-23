// Copyright (c) 2015, LAAS-CNRS
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

#ifndef HPP_CONSTRAINTS_CONVEX_HULL_HH
# define HPP_CONSTRAINTS_CONVEX_HULL_HH

# include <vector>
# include <hpp/constraints/differentiable-function.hh>
# include <hpp/fcl/math/transform.h>
# include <hpp/fcl/shape/geometric_shapes.h>
# include <hpp/constraints/tools.hh>

# include "hpp/constraints/fwd.hh"

namespace hpp {
  namespace constraints {
    /// Return the closest point to point P, on a segment line \f$ A + t*v, t \in [0,1] \f$.
    /// \param P PA where P is the point to
    /// \param A the origin the segment line
    /// \param v vector presenting the segment line
    /// \param[out] B the closest point
    inline void closestPointToSegment (const fcl::Vec3f& P,
        const fcl::Vec3f& A, const fcl::Vec3f& v, fcl::Vec3f& B) {
      fcl::Vec3f w = A - P;
      value_type c1, c2;
      c1 = v.dot (w);
      c2 = v.dot(v);
      if (c1 <= 0)       B = A;
      else if (c2 <= c1) B = A + v;
      else               B = A + c1 / c2 * v;
    }

    class HPP_CONSTRAINTS_DLLAPI ConvexShape
    {
      public:
        /// Represent a convex hull
        /// \param pts a sequence of points lying in a plane. The convex hull is
        ///        obtained by connecting consecutive points (in a circular way)
        ///
        /// \note There is no convexity check yet. The order is important:
        ///       The normal is parallel to (pts[1] - pts[0]).cross (pts[2] - pts[1])
        ///       The normal to the segment in the plane are directed outward.
        ///             (pts[i+1] - pts[i]).cross (normalToConvexShape)
        ConvexShape (const std::vector <vector3_t>& pts, JointPtr_t joint = NULL):
          Pts_ (pts), joint_ (joint)
        {
          init ();
        }

        ConvexShape (const fcl::TriangleP& t, const JointPtr_t& joint = NULL):
          Pts_ (triangleToPoints (t)), joint_ (joint)
        {
          init ();
        }

        // Copy constructor
        ConvexShape (const ConvexShape& t) :
          Pts_ (t.Pts_), joint_ (t.joint_)
        {
          init ();
        }

        void updateToCurrentTransform () const
        {
          if (joint_ != NULL) recompute (joint_->currentTransformation());
        }

        /// Intersection with a line defined by a point and a vector.
        /// updateToCurrentTransform() should be called before.
        inline fcl::Vec3f intersection (const fcl::Vec3f& A, const fcl::Vec3f& u) const {
          assert (std::abs (n_.dot (u)) > 1e-8);
          return A + u * (n_.dot(c_ - A)) / n_.dot (u);
        }
        inline fcl::Vec3f intersectionLocal (const fcl::Vec3f& A, const fcl::Vec3f& u) const {
          assert (std::abs (N_.dot (u)) > 1e-8);
          return A + u * (N_.dot(C_ - A)) / N_.dot (u);
        }

        /// Check whether the intersection of the line defined by A and u
        /// onto the plane containing the triangle is inside the triangle.
        inline bool isInside (const fcl::Vec3f& A, const fcl::Vec3f& u) const {
          assert (shapeDimension_ > 2);
          return isInside (intersection (A, u));
        }
        inline bool isInside (const fcl::Vec3f& Ap) const {
          assert (shapeDimension_ > 2);
          if (joint_ == NULL) return isInsideLocal (Ap);
          fcl::Transform3f M = joint_->currentTransformation ();
          vector3_t Ap_loc = M.inverse ().transform(Ap);
          return isInsideLocal (Ap_loc);
        }
        /// As isInside but consider A as expressed in joint frame.
        inline bool isInsideLocal (const fcl::Vec3f& Ap) const {
          assert (shapeDimension_ > 2);
          for (std::size_t i = 0; i < shapeDimension_; ++i) {
            if (Ns_[i].dot (Ap-Pts_[i]) > 0) return false;
          }
          return true;
        }

        /// Return the shortest distance from a point to the hull
        /// A negative value means the point is inside the hull
        /// \param A a point already in the plane containing the convex hull
        inline value_type distance (const fcl::Vec3f& A) const {
          assert (shapeDimension_ > 1);
          const value_type inf = std::numeric_limits<value_type>::infinity();
          value_type minPosDist = inf, maxNegDist = - inf;
          bool outside = false;
          for (std::size_t i = 0; i < shapeDimension_; ++i) {
            value_type d = dist (A - Pts_[i], Us_[i], Ns_[i]);
            if (d > 0) {
              outside = true;
              if (d < minPosDist) minPosDist = d;
            }
            if (d < 0 && d > maxNegDist) maxNegDist = d;
          }
          if (outside) return minPosDist;
          return maxNegDist;
        }

        inline const fcl::Vec3f& planeXaxis () const {
          assert (shapeDimension_ > 2);
          return   n0_;
        }
        inline const fcl::Vec3f& planeYaxis () const {
          assert (shapeDimension_ > 2);
          return nxn0_;
        }
        inline const fcl::Vec3f& normal () const {
          assert (shapeDimension_ > 2);
          return n_;
        }
        inline const fcl::Vec3f& center () const { return c_; }
        /// Transform from world frame coordinate to local frame coordinate
        inline const fcl::Transform3f& inversePosition () const { return M_; }

        /// The position in the joint frame and the joint
        std::vector <vector3_t> Pts_;
        size_t shapeDimension_;
        vector3_t C_, N_;
        std::vector <vector3_t> Ns_, Us_;
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
          c2 = v.dot (v);
          if (c2 <= c1)
            return (u.dot (w) > 0)?((w-v).norm()):(-(w-v).norm());
          return u.dot (w);
        }

        static std::vector <vector3_t> triangleToPoints (const fcl::TriangleP& t) {
          std::vector <vector3_t> ret (3);
          ret[0] = t.a;
          ret[1] = t.b;
          ret[2] = t.c;
          return ret;
        }

        void init ()
        {
          shapeDimension_ = Pts_.size ();

          switch (shapeDimension_) {
            case 0:
              throw std::logic_error ("Cannot represent an empty shape.");
              break;
            case 1:
              C_ = Pts_[0];
              // The transformation will be (N_, Ns_[0], Us_[0])
              // Fill vectors so as to be consistent
              N_ = vector3_t(1,0,0);
              Ns_.push_back (vector3_t(0,1,0));
              Us_.push_back (vector3_t(0,0,1));
              break;
            case 2:
              C_ = (Pts_[0] + Pts_[1])/2;
              // The transformation will be (N_, Ns_[0], Us_[0])
              // Fill vectors so as to be consistent
              Us_.push_back (Pts_[1] - Pts_[0]);
              Us_[0].normalize ();
              if (Us_[0][0] != 0) N_ = vector3_t(-Us_[0][1],Us_[0][0],0);
              else                N_ = vector3_t(0,-Us_[0][2],Us_[0][1]);
              N_.normalize ();
              Ns_.push_back (Us_[0].cross (N_));
              Ns_[0].normalize (); // Should be unnecessary
              break;
            default:
              C_.setZero ();
              for (std::size_t i = 0; i < shapeDimension_; ++i)
                C_ += Pts_[i];
              C_ /= Pts_.size();
              N_ = (Pts_[1] - Pts_[0]).cross (Pts_[2] - Pts_[1]);
              assert (!N_.isZero ());
              N_.normalize ();

              Us_.resize (Pts_.size());
              Ns_.resize (Pts_.size());
              for (std::size_t i = 0; i < shapeDimension_; ++i) {
                Us_[i] = Pts_[(i+1)%shapeDimension_] - Pts_[i];
                Us_[i].normalize ();
                Ns_[i] = Us_[i].cross (N_);
                Ns_[i].normalize ();
              }
              for (std::size_t i = 0; i < shapeDimension_; ++i) {
                assert (Us_[(i+1)%shapeDimension_].dot (Ns_[i]) < 0 &&
                    "The sequence does not define a convex surface");
              }
              break;
          }

          if (joint_ == NULL) recompute (Transform3f ());
          else                recompute (joint_->currentTransformation ());
        }

        void recompute (const Transform3f& M) const
        {
          c_ = M.transform (C_);
          n_ = M.getRotation () * (N_);
          n0_ = M.getRotation () * (Ns_[0]);
          nxn0_ = M.getRotation () * (Us_[0]);
          fcl::Matrix3f R (n_, n0_, nxn0_);
          M_ = fcl::Transform3f (R, - (R * c_) );
          // M_.setTranslation (- (M_.getRotation () * c_));
          for (size_t i = 0; i < 3; i++) assert (M_.getRotation () (0, i) == n_[i]);
          for (size_t i = 0; i < 3; i++) assert (M_.getRotation () (1, i) == n0_[i]);
          for (size_t i = 0; i < 3; i++) assert (M_.getRotation () (2, i) == nxn0_[i]);
        }

        /// The positions and vectors in the global frame
        mutable fcl::Vec3f n_, c_, n0_, nxn0_;
        mutable fcl::Transform3f M_;
    };
  } // namespace constraints
} // namespace hpp

#endif //  HPP_CONSTRAINTS_CONVEX_HULL_HH
