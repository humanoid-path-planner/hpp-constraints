// Copyright (c) 2015, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
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

#define BOOST_TEST_MODULE ConvexShape

#include "hpp/constraints/convex-shape.hh"

#include <boost/test/included/unit_test.hpp>
#include <pinocchio/fwd.hpp>

using hpp::constraints::ConvexShape;
using hpp::constraints::ConvexShapeData;
using hpp::constraints::value_type;
using hpp::constraints::vector3_t;

BOOST_AUTO_TEST_CASE(triangle) {
  vector3_t p0(0, 0, 0), p1(2, 0, 0), p2(0, 2, 0);
  std::vector<vector3_t> pts;
  pts.push_back(p0);
  pts.push_back(p1);
  pts.push_back(p2);
  ConvexShape t(pts);
  ConvexShapeData d;
  d.updateToCurrentTransform(t);

  std::vector<vector3_t> ptsIn, ptsOut;
  ptsIn.push_back(p0);
  ptsIn.push_back(p1);
  ptsIn.push_back(p2);
  ptsIn.push_back(vector3_t(1, 1, 0));
  ptsIn.push_back(vector3_t(0.5, 0.5, 0));

  ptsOut.push_back(vector3_t(-1, -1, 0));
  ptsOut.push_back(vector3_t(1, -1, 0));
  ptsOut.push_back(vector3_t(3, -0.2, 0));
  ptsOut.push_back(vector3_t(2, 2, 0));
  ptsOut.push_back(vector3_t(-0.2, 3, 0));
  ptsOut.push_back(vector3_t(-1, 1, 0));

  for (size_t i = 0; i < ptsIn.size(); i++) {
    BOOST_CHECK_EQUAL(t.isInsideLocal(ptsIn[i]), d.isInside(t, ptsIn[i]));
    BOOST_CHECK_EQUAL(t.distanceLocal(ptsIn[i]), d.distance(t, ptsIn[i]));

    BOOST_CHECK_MESSAGE(
        d.isInside(t, ptsIn[i]),
        "Check point inside failed for ptsIn[" << i << "]=" << ptsIn[i]);
    BOOST_CHECK_MESSAGE(d.distance(t, ptsIn[i]) <= 0,
                        "Wrong point to triangle distance for ptsIn["
                            << i << "]=" << ptsIn[i]
                            << ". Distance returned is "
                            << d.distance(t, ptsIn[i]));
  }
  for (size_t i = 0; i < ptsOut.size(); i++) {
    BOOST_CHECK_MESSAGE(
        !d.isInside(t, ptsOut[i]),
        "Check point outside failed for ptsOut[" << i << "]=" << ptsOut[i]);
    BOOST_CHECK_MESSAGE(d.distance(t, ptsOut[i]) >= 0,
                        "Wrong point to triangle distance for ptsOut["
                            << i << "]=" << ptsOut[i]
                            << ". Distance returned is "
                            << d.distance(t, ptsOut[i]));
  }
}

void checkDistance(const ConvexShape& t, const vector3_t& p,
                   const value_type& expected) {
  ConvexShapeData d;
  d.updateToCurrentTransform(t);

  BOOST_CHECK_MESSAGE(std::abs(d.distance(t, p) - expected) < 1e-5,
                      "Point " << p
                               << ".\nDistance computed: " << d.distance(t, p)
                               << "\nDistance expected: " << expected);
}

BOOST_AUTO_TEST_CASE(distance) {
  vector3_t p0(0, 0, 0), p1(2, 0, 0), p2(2, 2, 0), p3(0, 2, 0);
  std::vector<vector3_t> pts;
  pts.push_back(p0);
  pts.push_back(p1);
  pts.push_back(p2);
  pts.push_back(p3);
  ConvexShape t(pts);

  checkDistance(t, vector3_t(0, 3, 0), 1);
  checkDistance(t, vector3_t(3, 0, 0), 1);
  checkDistance(t, vector3_t(1, 1, 0), -1);
  checkDistance(t, vector3_t(0, 1, 0), 0);
}
