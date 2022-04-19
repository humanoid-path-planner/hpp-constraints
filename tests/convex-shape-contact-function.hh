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

#include <hpp/constraints/convex-shape-contact.hh>

hpp::constraints::ConvexShapeContactPtr_t createConvexShapeContact_punctual(
    hpp::pinocchio::DevicePtr_t d, hpp::pinocchio::JointPtr_t j,
    std::string name) {
  using namespace hpp::constraints;
  /** Floor = penta + square
   *     +
   *    / \
   *   /   +-----+
   *  +    |     |
   *   \   +-----+
   *    \ /
   *     +
   *
   *  Object = point
   **/
  std::vector<vector3_t> square(4);
  square[0] = vector3_t(5, 5, 0);
  square[1] = vector3_t(5, -5, 0);
  square[2] = vector3_t(0, -5, 0);
  square[3] = vector3_t(0, 5, 0);

  std::vector<vector3_t> penta(5);
  penta[0] = vector3_t(0, 5, 0);
  penta[1] = vector3_t(0, -5, 0);
  penta[2] = vector3_t(-2, -6, 0);
  penta[3] = vector3_t(-5, 0, 0);
  penta[4] = vector3_t(-2, 6, 0);

  std::vector<vector3_t> point(1, vector3_t(0, 0, 0));

  JointAndShapes_t os;
  os.push_back(JointAndShape_t(j, point));
  JointAndShapes_t fs;
  fs.push_back(JointAndShape_t(JointPtr_t(), square));
  fs.push_back(JointAndShape_t(JointPtr_t(), penta));
  ConvexShapeContactPtr_t fptr = ConvexShapeContact::create(name, d, fs, os);
  return fptr;
}

hpp::constraints::ConvexShapeContactPtr_t createConvexShapeContact_convex(
    hpp::pinocchio::DevicePtr_t d, hpp::pinocchio::JointPtr_t j,
    std::string name) {
  using namespace hpp::constraints;
  /** Floor = penta + square
   *     +
   *    / \
   *   /   +-----+
   *  +    |     |
   *   \   +-----+
   *    \ /
   *     +
   *
   *  Object = trapezium
   *  +--+
   *  |   \
   *  +----+
   **/
  std::vector<vector3_t> square(4);
  square[0] = vector3_t(5, 5, 0);
  square[1] = vector3_t(5, -5, 0);
  square[2] = vector3_t(0, -5, 0);
  square[3] = vector3_t(0, 5, 0);

  std::vector<vector3_t> penta(5);
  penta[0] = vector3_t(0, 5, 0);
  penta[1] = vector3_t(0, -5, 0);
  penta[2] = vector3_t(-2, -6, 0);
  penta[3] = vector3_t(-5, 0, 0);
  penta[4] = vector3_t(-2, 6, 0);

  std::vector<vector3_t> trapeze(4);
  trapeze[0] = vector3_t(-0.1, 0.1, 0);
  trapeze[1] = vector3_t(0.1, 0.1, 0);
  trapeze[2] = vector3_t(0.2, -0.1, 0);
  trapeze[3] = vector3_t(-0.1, -0.1, 0);
  ConvexShape tCs(trapeze, j);

  JointAndShapes_t os;
  os.push_back(JointAndShape_t(j, trapeze));
  JointAndShapes_t fs;
  fs.push_back(JointAndShape_t(JointPtr_t(), square));
  fs.push_back(JointAndShape_t(JointPtr_t(), penta));
  ConvexShapeContactPtr_t fptr = ConvexShapeContact::create(name, d, fs, os);
  return fptr;
}

hpp::constraints::ConvexShapeContactPtr_t createConvexShapeContact_triangles(
    hpp::pinocchio::DevicePtr_t d, hpp::pinocchio::JointPtr_t j,
    std::string name) {
  using namespace hpp::constraints;
  /** Floor = penta + square
   *     +
   *    / \
   *   /   +-----+
   *  +    |     |
   *   \   +-----+
   *    \ /
   *     +
   *
   *  Object = point
   *  +--+
   *  |   \
   *  +----+
   **/
  std::vector<vector3_t> square(4);
  square[0] = vector3_t(5, 5, 0);
  square[1] = vector3_t(5, -5, 0);
  square[2] = vector3_t(0, -5, 0);
  square[3] = vector3_t(0, 5, 0);

  std::vector<vector3_t> penta(5);
  penta[0] = vector3_t(0, 5, 0);
  penta[1] = vector3_t(0, -5, 0);
  penta[2] = vector3_t(-2, -6, 0);
  penta[3] = vector3_t(-5, 0, 0);
  penta[4] = vector3_t(-2, 6, 0);

  std::vector<vector3_t> trapeze(4);
  trapeze[0] = vector3_t(-0.1, 0.1, 0);
  trapeze[1] = vector3_t(0.1, 0.1, 0);
  trapeze[2] = vector3_t(0.2, -0.1, 0);
  trapeze[3] = vector3_t(-0.1, -0.1, 0);

  JointAndShapes_t os;
  os.push_back(JointAndShape_t(j, trapeze));
  JointAndShapes_t fs;
  fs.push_back(JointAndShape_t(JointPtr_t(), square));
  fs.push_back(JointAndShape_t(JointPtr_t(), penta));
  ConvexShapeContactPtr_t fptr = ConvexShapeContact::create(name, d, fs, os);
  return fptr;
}
