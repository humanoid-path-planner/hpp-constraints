// Copyright (c) 2014, LAAS-CNRS
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

#define BOOST_TEST_MODULE ConvexShapes
#include <boost/test/included/unit_test.hpp>

#include "hpp/constraints/convex-shape.hh"

using hpp::constraints::ConvexShape;

BOOST_AUTO_TEST_CASE(triangleAsConvexShape) {
  coal::Vec3f p0(0, 0, 0), p1(2, 0, 0), p2(0, 2, 0);
  ConvexShape t(p0, p1, p2);

  std::vector<coal::Vec3f> ptsIn, ptsOut;
  ptsIn.push_back(p0);
  ptsIn.push_back(p1);
  ptsIn.push_back(p2);
  ptsIn.push_back(coal::Vec3f(1, 1, 0));
  ptsIn.push_back(coal::Vec3f(0.5, 0.5, 0));

  ptsOut.push_back(coal::Vec3f(-1, -1, 0));
  ptsOut.push_back(coal::Vec3f(1, -1, 0));
  ptsOut.push_back(coal::Vec3f(3, -0.2, 0));
  ptsOut.push_back(coal::Vec3f(2, 2, 0));
  ptsOut.push_back(coal::Vec3f(-0.2, 3, 0));
  ptsOut.push_back(coal::Vec3f(-1, 1, 0));

  for (size_t i = 0; i < ptsIn.size(); i++) {
    BOOST_CHECK_MESSAGE(
        t.isInside(ptsIn[i]),
        "Check point inside failed for ptsIn[" << i << "]=" << ptsIn[i]);
    BOOST_CHECK_MESSAGE(t.distance(ptsIn[i]) <= 0,
                        "Wrong point to triangle distance for ptsIn["
                            << i << "]=" << ptsIn[i]
                            << ". Distance returned is "
                            << t.distance(ptsIn[i]));
  }
  for (size_t i = 0; i < ptsOut.size(); i++) {
    BOOST_CHECK_MESSAGE(
        !t.isInside(ptsOut[i]),
        "Check point inside failed for ptsOut[" << i << "]=" << ptsOut[i]);
    BOOST_CHECK_MESSAGE(t.distance(ptsOut[i]) >= 0,
                        "Wrong point to triangle distance for ptsOut["
                            << i << "]=" << ptsOut[i]
                            << ". Distance returned is "
                            << t.distance(ptsOut[i]));
  }
}
