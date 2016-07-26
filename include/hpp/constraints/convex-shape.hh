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

#ifndef HPP_CONSTRAINTS_CONVEX_SHAPE_HH
# define HPP_CONSTRAINTS_CONVEX_SHAPE_HH

# include <vector>

# include <hpp/fcl/math/transform.h>
# include <hpp/fcl/shape/geometric_shapes.h>

# include <hpp/pinocchio/joint.hh>

# include <hpp/constraints/differentiable-function.hh>

namespace hpp {
  namespace constraints {
    /// Return the closest point to point P, on a segment line \f$ A + t*v, t \in [0,1] \f$.
    /// \param P PA where P is the point to
    /// \param A the origin the segment line
    /// \param v vector presenting the segment line
    /// \param[out] B the closest point
    inline void closestPointToSegment (const vector3_t& P,
        const vector3_t& A, const vector3_t& v, vector3_t& B) {
      vector3_t w = A - P;
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
        /// Represent a convex shape
        /// \param pts a sequence of points lying in a plane. The convex shape is
        ///        obtained by connecting consecutive points (in a circular way)
        ///
        /// \note There is no convexity check yet. The order is important:
        ///       The normal is parallel to (pts[1] - pts[0]).cross (pts[2] - pts[1])
        ///       The normal to the segment in the plane are directed outward.
        ///             (pts[i+1] - pts[i]).cross (normalToConvexShape)
        ConvexShape (const std::vector <vector3_t>& pts, JointPtr_t joint = JointPtr_t()):
          Pts_ (pts), joint_ (joint)
        {
          init ();
        }

        ConvexShape (const fcl::TriangleP& t, const JointPtr_t& joint = JointPtr_t()):
          Pts_ (triangleToPoints (t)), joint_ (joint)
        {
          init ();
        }

        /// This constructor is required for compatibility with deprecated
        /// Triangle constructor.
        ConvexShape (const vector3_t& p0, const vector3_t& p1,
            const vector3_t& p2, const JointPtr_t& joint = JointPtr_t()):
          Pts_ (points(p0,p1,p2)), joint_ (joint)
        {
          init ();
        }

        // Copy constructor
        ConvexShape (const ConvexShape& t) :
          Pts_ (t.Pts_), joint_ (t.joint_)
        {
          init ();
        }

        void reverse () {
          std::reverse (Pts_.begin (), Pts_.end());
          init ();
        }

        void updateToCurrentTransform () const
        {
          if (joint_ != NULL) recompute (joint_->currentTransformation());
        }

        /// Intersection with a line defined by a point and a vector.
        /// updateToCurrentTransform() should be called before.
        inline vector3_t intersection (const vector3_t& A, const vector3_t& u) const {
          assert (std::abs (n_.dot (u)) > 1e-8);
          return A + u * (n_.dot(c_ - A)) / n_.dot (u);
        }
        inline vector3_t intersectionLocal (const vector3_t& A, const vector3_t& u) const {
          assert (std::abs (N_.dot (u)) > 1e-8);
          return A + u * (N_.dot(C_ - A)) / N_.dot (u);
        }

        /// Check whether the intersection of the line defined by A and u
        /// onto the plane containing the triangle is inside the triangle.
        inline bool isInside (const vector3_t& A, const vector3_t& u) const {
          assert (shapeDimension_ > 2);
          return isInside (intersection (A, u));
        }
        inline bool isInside (const vector3_t& Ap) const {
          assert (shapeDimension_ > 2);
          if (joint_ == NULL) return isInsideLocal (Ap);
          const Transform3f& M = joint_->currentTransformation ();
          vector3_t Ap_loc = M.actInv(Ap);
          return isInsideLocal (Ap_loc);
        }
        /// As isInside but consider A as expressed in joint frame.
        inline bool isInsideLocal (const vector3_t& Ap) const {
          assert (shapeDimension_ > 2);
          for (std::size_t i = 0; i < shapeDimension_; ++i) {
            if (Ns_[i].dot (Ap-Pts_[i]) > 0) return false;
          }
          return true;
        }

        /// Return the shortest distance from a point to the shape
        /// A negative value means the point is inside the shape
        /// \param A a point already in the plane containing the convex shape,
        ///        and expressed in the global frame.
        inline value_type distance (vector3_t a) const {
          assert (shapeDimension_ > 1);
          if (joint_!=NULL) a = joint_->currentTransformation ().actInv(a);
          const value_type inf = std::numeric_limits<value_type>::infinity();
          value_type minPosDist = inf, maxNegDist = - inf;
          bool outside = false;
          for (std::size_t i = 0; i < shapeDimension_; ++i) {
            value_type d = dist (a - Pts_[i], Ls_[i], Us_[i], Ns_[i]);
            if (d > 0) {
              outside = true;
              if (d < minPosDist) minPosDist = d;
            }
            if (d <= 0 && d > maxNegDist) maxNegDist = d;
          }
          if (outside) return minPosDist;
          return maxNegDist;
        }

        /// Return the X axis of the plane in the joint frame
        inline const vector3_t& planeXaxis () const {
          assert (shapeDimension_ > 2);
          return Ns_[0];
        }
        /// Return the Y axis of the plane in the joint frame
        /// The Y axis is aligned with \f$ Pts_[1] - Pts_[0] \f$
        inline const vector3_t& planeYaxis () const {
          assert (shapeDimension_ > 2);
          return Us_[0];
        }

        /// Return the normal in world frame.
        inline const vector3_t& normal () const {
          assert (shapeDimension_ > 2);
          return n_;
        }
        /// Return the center in world frame.
        inline const vector3_t& center () const { return c_; }

        /// Transform of the shape in the joint frame
        inline const Transform3f& positionInJoint () const { return MinJoint_; }
        // \param yaxis vector in world frame to which we should try to align
        inline void computeAlignedPosition (vector3_t yaxis) const {
          assert (shapeDimension_ > 2);
          // Project vector onto the plane
          if (joint_!=NULL) yaxis = joint_->currentTransformation ().actInv(yaxis);
          vector3_t yproj = yaxis - yaxis.dot (N_) * N_;
          if (yproj.isZero ()) M_ = MinJoint_;
          else {
            M_.translation() = C_;
            M_.rotation().col(0) = N_;
            M_.rotation().col(1) = yaxis;
            M_.rotation().col(2) = N_.cross (yaxis);
          }
        }
        inline const Transform3f& alignedPositionInJoint () const { return M_; }

        /// The points in the joint frame. It is constant.
        std::vector <vector3_t> Pts_;
        size_t shapeDimension_;
        /// the center in the joint frame. It is constant.
        vector3_t C_;
        /// the normal to the shape in the joint frame. It is constant.
        vector3_t N_;
        /// Ns_ and Us_ are unit vector, in the plane containing the shape,
        /// expressed in the joint frame.
        /// Ns_[i] is normal to edge i, pointing inside.
        /// Ns_[i] is a vector director of edge i.
        std::vector <vector3_t> Ns_, Us_;
        vector_t Ls_;
        Transform3f MinJoint_;
        JointPtr_t joint_;

      private:
        /// Return the distance between the point A and the segment
        /// [P, c2*v] oriented by u.
        /// w = PA.
        inline value_type dist (const vector3_t& w, const value_type& c2, const vector3_t& v, const vector3_t& u) const {
          value_type c1;
          c1 = v.dot (w);
          if (c1 <= 0)
            return (u.dot (w) > 0)?(w.norm()):(- w.norm());
          if (c2 <= c1)
            // TODO: (w - c2 * v).norm() == sqrt((u.dot(w)**2 + (c1 - c2)**2)
            // second should be cheaper.
            return (u.dot (w) > 0)?((w-c2*v).norm()):(-(w-c2*v).norm());
          return u.dot (w);
        }

        static std::vector <vector3_t> triangleToPoints (const fcl::TriangleP& t) {
          // TODO
          // return points (t.a, t.b, t.c);
          std::vector <vector3_t> ret (3);
          ret[0] = t.a;
          ret[1] = t.b;
          ret[2] = t.c;
          return ret;
        }
        static std::vector <vector3_t> points (const vector3_t& p0,
            const vector3_t& p1, const vector3_t& p2) {
          std::vector <vector3_t> ret (3);
          ret[0] = p0; ret[1] = p1; ret[2] = p2;
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
              Ls_ = vector_t(1);
              C_ = (Pts_[0] + Pts_[1])/2;
              // The transformation will be (N_, Ns_[0], Us_[0])
              // Fill vectors so as to be consistent
              Us_.push_back (Pts_[1] - Pts_[0]);
              Ls_[0] = Us_[0].norm();
              Us_[0].normalize ();
              if (Us_[0][0] != 0) N_ = vector3_t(-Us_[0][1],Us_[0][0],0);
              else                N_ = vector3_t(0,-Us_[0][2],Us_[0][1]);
              N_.normalize ();
              Ns_.push_back (Us_[0].cross (N_));
              Ns_[0].normalize (); // Should be unnecessary
              break;
            default:
              Ls_ = vector_t(shapeDimension_);
              C_.setZero ();
              for (std::size_t i = 0; i < shapeDimension_; ++i)
                C_ += Pts_[i];
              // TODO This is very ugly. Why Eigen does not have the operator/=(int) ...
              C_ /= (value_type)Pts_.size();
              N_ = (Pts_[1] - Pts_[0]).cross (Pts_[2] - Pts_[1]);
              assert (!N_.isZero ());
              N_.normalize ();

              Us_.resize (Pts_.size());
              Ns_.resize (Pts_.size());
              for (std::size_t i = 0; i < shapeDimension_; ++i) {
                Us_[i] = Pts_[(i+1)%shapeDimension_] - Pts_[i];
                Ls_[i] = Us_[i].norm();
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

          MinJoint_.translation() = C_;
          MinJoint_.rotation().col(0) = N_;
          MinJoint_.rotation().col(1) = Ns_[0];
          MinJoint_.rotation().col(2) = Us_[0];

          if (joint_ == NULL) recompute (Transform3f::Identity());
          else                recompute (joint_->currentTransformation ());
        }

        void recompute (const Transform3f& M) const
        {
          c_ = M.act (C_);
          n_ = M.rotation () * N_;
        }

        /// The positions and vectors in the global frame
        mutable vector3_t n_, c_;
        mutable Transform3f M_;
    };
  } // namespace constraints
} // namespace hpp

#endif //  HPP_CONSTRAINTS_CONVEX_SHAPE_HH
