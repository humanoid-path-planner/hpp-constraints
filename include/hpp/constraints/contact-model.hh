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

#ifndef HPP_CONSTRAINTS_CONTACT_HH
#define HPP_CONSTRAINTS_CONTACT_HH

#include <boost/assign/list_of.hpp>

#include "hpp/constraints/config.hh"
#include "hpp/constraints/fwd.hh"

#include "hpp/constraints/convex-shape.hh"

namespace hpp {
  namespace constraints {
    struct ForceData {
      JointPtr_t j;
      std::vector<vector3_t> p;
      vector3_t n;
      //vector3_t tX, tY;
    };

    class PointSurface;
    class LineSurface;
    class FlatSurface;

    // Square, triangle, segment or point
    class HPP_CONSTRAINTS_DLLAPI AbstractSurface {
      public:
        // typedef Eigen::Matrix <value_type, 3, 1> Point_t;
        typedef vector3_t Point_t;
        typedef std::vector <Point_t> Points_t;
        typedef std::vector <ForceData> ForceDatas_t;

        AbstractSurface (Points_t pts, JointPtr_t j) : 
          cs_ (pts, j), joint_ (j) {}

        virtual Point_t center () const = 0;

        // Returns the normal direction between two surface
        Point_t normal (const AbstractSurface* other) const;
        virtual Point_t normal (const PointSurface* other) const = 0;
        virtual Point_t normal (const LineSurface* other) const = 0;
        virtual Point_t normal (const FlatSurface* other) const = 0;

        // Axis for the contact forces
        // virtual Point_t tangentX () const = 0;
        // virtual Point_t tangentY () const = 0;

        // Project point p onto the surface.
        virtual Point_t project (const Point_t& p) const = 0;

        // Returns a list of ForceData representing the interaction between the
        // two surfaces.
        virtual ForceDatas_t forcePoints (const AbstractSurface* other) const = 0;
        // virtual Points_t forcePointsInCurrentPosition ()

        // Returns the friction coefficient between the two surface
        virtual value_type frictionCoeff (AbstractSurface* /*other*/) const
        { return 0.7; }

        // Transform point to global frame coordinates
        Point_t pointToGlobalFrame  (const Point_t& p) const;
        Point_t vectorToGlobalFrame (const Point_t& p) const;
        Point_t pointToLocalFrame   (const Point_t& p) const;
        Point_t vectorToLocalFrame  (const Point_t& p) const;

        JointPtr_t joint () const {
          return joint_;
        }

        // Get the underlying convex shape
        const ConvexShape& convexShape () const;

      protected:
        ConvexShape cs_;
        JointPtr_t joint_;
    };

    // Has to be seen as a sphere of radius 0
    class HPP_CONSTRAINTS_DLLAPI PointSurface : public AbstractSurface {
      public:
        PointSurface (Point_t point, JointPtr_t joint = NULL) :
          AbstractSurface (Points_t (1, point), joint), p_ (point)
        {}

        Point_t center () const { return p_; }

        Point_t project (const Point_t& /*p*/) const { return p_; }

        Point_t normal (const PointSurface* other) const;
        Point_t normal (const LineSurface* other) const;
        Point_t normal (const FlatSurface* other) const;

        // Returns a list a points where the forces are applied.
        ForceDatas_t forcePoints (const AbstractSurface* other) const;

      private:
        Point_t p_, n_;
    };

    // Has to be seen as a cylinder of radius 0
    class HPP_CONSTRAINTS_DLLAPI LineSurface : public AbstractSurface {
      public:
        LineSurface (Point_t p0, Point_t p1, JointPtr_t joint = NULL) :
          AbstractSurface (boost::assign::list_of (p0)(p1), joint),
          p0_ (p0), p1_ (p1)
        {
          assert ((p0_ - p1_).isZero ());
        }

        Point_t center () const { return 0.5 * (p0_ + p1_); }

        Point_t project (const Point_t& p) const
        {
          const fcl::Vec3f u = p1_ - p0_;
          return p0_ + u * (p - p0_).dot(u) / u.dot (u);
        }

        Point_t normal (const PointSurface* other) const;
        Point_t normal (const LineSurface* other) const;
        Point_t normal (const FlatSurface* other) const;

        // Returns a list a points where the forces are applied.
        ForceDatas_t forcePoints (const AbstractSurface* other) const;

      private:
        Point_t p0_, p1_;
        JointPtr_t j_;
    };

    class HPP_CONSTRAINTS_DLLAPI FlatSurface : public AbstractSurface {
      public:
        FlatSurface (Points_t pts, JointPtr_t joint = NULL) :
          AbstractSurface (pts, joint)
        {}

        Point_t normal (const PointSurface* /*other*/) const { return normal(); }
        Point_t normal (const LineSurface*  /*other*/) const { return normal(); }
        Point_t normal (const FlatSurface*  /*other*/) const { return normal(); }

        Point_t center () const { return cs_.C_; }
        Point_t normal () const { return cs_.N_; }

        Point_t tangentX () const { return cs_.Ns_[0]; }
        Point_t tangentY () const { return cs_.Us_[0]; }

        Point_t project (const Point_t& p) const {
          return cs_.intersectionLocal (p, cs_.N_);
        };

        // Returns a list a points where the forces are applied.
        ForceDatas_t forcePoints (const AbstractSurface* other) const;
    };

    // Returns a constraint ensuring a null velocity between the surface
    DifferentiableFunctionPtr_t nullVelocityConstraint
      (const FlatSurface& s1, const FlatSurface& s2, const DevicePtr_t device);
    DifferentiableFunctionPtr_t nullVelocityConstraint
      (const PointSurface& s1, const FlatSurface& s2, const DevicePtr_t device);

    // TODO: this should be handled in nullVelocityConstraint with a boolean
    // noSlide for projection efficiency.
    DifferentiableFunctionPtr_t noSlideConstraint
      (const FlatSurface& s1, const FlatSurface& s2, const DevicePtr_t device);
    DifferentiableFunctionPtr_t noSlideConstraint
      (const PointSurface& s1, const FlatSurface& s2, const DevicePtr_t device);

    // Returns a constraint that ensure the two surfaces are in contact
    // nullVelocityConstraint + constraint along XY axis of plane s2
    DifferentiableFunctionPtr_t contactConstraint
      (const FlatSurface& s1, const FlatSurface& s2, const DevicePtr_t device);
    DifferentiableFunctionPtr_t contactConstraint
      (const PointSurface& s1, const FlatSurface& s2, const DevicePtr_t device);
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_CONTACT_HH
