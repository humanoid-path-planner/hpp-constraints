// Copyright (c) 2015, Joseph Mirabel
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

#define BOOST_TEST_MODULE ConvexShape
#include <boost/test/included/unit_test.hpp>

#include "hpp/constraints/convex-shape.hh"

using hpp::constraints::ConvexShape;
using hpp::constraints::value_type;

BOOST_AUTO_TEST_CASE (triangle)
{
  fcl::Vec3f p0 (0,0,0),
             p1 (2,0,0),
             p2 (0,2,0);
  std::vector <fcl::Vec3f> pts;
  pts.push_back (p0);
  pts.push_back (p1);
  pts.push_back (p2);
  ConvexShape t (pts);

  std::vector <fcl::Vec3f> ptsIn, ptsOut;
  ptsIn.push_back (p0);
  ptsIn.push_back (p1);
  ptsIn.push_back (p2);
  ptsIn.push_back (fcl::Vec3f (1,1,0));
  ptsIn.push_back (fcl::Vec3f (0.5,0.5,0));

  ptsOut.push_back (fcl::Vec3f (-1,-1,0));
  ptsOut.push_back (fcl::Vec3f (1,-1,0));
  ptsOut.push_back (fcl::Vec3f (3,-0.2,0));
  ptsOut.push_back (fcl::Vec3f ( 2,2,0));
  ptsOut.push_back (fcl::Vec3f (-0.2,3,0));
  ptsOut.push_back (fcl::Vec3f (-1,1,0));

  for (size_t i = 0; i < ptsIn.size (); i++) {
    BOOST_CHECK_MESSAGE (t.isInside (ptsIn[i]),
        "Check point inside failed for ptsIn[" << i << "]=" << ptsIn[i]);
    BOOST_CHECK_MESSAGE (t.distance (ptsIn[i]) <= 0,
        "Wrong point to triangle distance for ptsIn[" << i << "]=" << ptsIn[i]
        << ". Distance returned is " << t.distance (ptsIn[i]));
  }
  for (size_t i = 0; i < ptsOut.size (); i++) {
    BOOST_CHECK_MESSAGE (!t.isInside (ptsOut[i]),
        "Check point outside failed for ptsOut[" << i << "]=" << ptsOut[i]);
    BOOST_CHECK_MESSAGE (t.distance (ptsOut[i]) >= 0,
        "Wrong point to triangle distance for ptsOut[" << i << "]=" << ptsOut[i]
        << ". Distance returned is " << t.distance (ptsOut[i]));
  }
}

void checkDistance(const ConvexShape& t, const fcl::Vec3f& p, const value_type& expected)
{
  BOOST_CHECK_MESSAGE(std::abs(t.distance(p) - expected) < 1e-5,
      "Point " << p << ".\nDistance computed: " << t.distance(p) << "\nDistance expected: " << expected
      );
}

BOOST_AUTO_TEST_CASE (distance)
{
  fcl::Vec3f p0 (0,0,0),
             p1 (2,0,0),
             p2 (2,2,0),
             p3 (0,2,0);
  std::vector <fcl::Vec3f> pts;
  pts.push_back (p0);
  pts.push_back (p1);
  pts.push_back (p2);
  pts.push_back (p3);
  ConvexShape t (pts);

  checkDistance(t, fcl::Vec3f(0, 3, 0),  1);
  checkDistance(t, fcl::Vec3f(3, 0, 0),  1);
  checkDistance(t, fcl::Vec3f(1, 1, 0), -1);
  checkDistance(t, fcl::Vec3f(0, 1, 0),  0);
}
