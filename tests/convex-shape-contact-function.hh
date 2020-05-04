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

#include <hpp/constraints/convex-shape-contact.hh>

hpp::constraints::ConvexShapeContactPtr_t createConvexShapeContact_punctual (hpp::pinocchio::DevicePtr_t d, hpp::pinocchio::JointPtr_t j, std::string name)
{
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
  std::vector <vector3_t> square(4);
  square[0] = vector3_t ( 5, 5,0); square[1] = vector3_t ( 5,-5,0);
  square[2] = vector3_t ( 0,-5,0); square[3] = vector3_t ( 0, 5,0);

  std::vector <vector3_t> penta(5);
  penta[0] = vector3_t ( 0, 5,0); penta[1] = vector3_t ( 0,-5,0);
  penta[2] = vector3_t (-2,-6,0); penta[3] = vector3_t (-5, 0,0);
  penta[4] = vector3_t (-2, 6,0);

  std::vector <vector3_t> point(1, vector3_t (0,0,0));

  JointAndShapes_t os; os.push_back(JointAndShape_t(j, point));
  JointAndShapes_t fs;
  fs.push_back(JointAndShape_t(JointPtr_t(),square));
  fs.push_back(JointAndShape_t(JointPtr_t(),penta));
  ConvexShapeContactPtr_t fptr = ConvexShapeContact::create (name, d, fs, os);
  return fptr;
}

hpp::constraints::ConvexShapeContactPtr_t createConvexShapeContact_convex (hpp::pinocchio::DevicePtr_t d, hpp::pinocchio::JointPtr_t j, std::string name)
{
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
  std::vector <vector3_t> square(4);
  square[0] = vector3_t ( 5, 5,0); square[1] = vector3_t ( 5,-5,0);
  square[2] = vector3_t ( 0,-5,0); square[3] = vector3_t ( 0, 5,0);

  std::vector <vector3_t> penta(5);
  penta[0] = vector3_t ( 0, 5,0); penta[1] = vector3_t ( 0,-5,0);
  penta[2] = vector3_t (-2,-6,0); penta[3] = vector3_t (-5, 0,0);
  penta[4] = vector3_t (-2, 6,0);

  std::vector <vector3_t> trapeze(4);
  trapeze[0] = vector3_t (-0.1, 0.1,0); trapeze[1] = vector3_t ( 0.1, 0.1,0);
  trapeze[2] = vector3_t ( 0.2,-0.1,0); trapeze[3] = vector3_t (-0.1,-0.1,0);
  ConvexShape tCs (trapeze, j);

  JointAndShapes_t os; os.push_back(JointAndShape_t(j, trapeze));
  JointAndShapes_t fs;
  fs.push_back(JointAndShape_t(JointPtr_t(),square));
  fs.push_back(JointAndShape_t(JointPtr_t(),penta));
  ConvexShapeContactPtr_t fptr = ConvexShapeContact::create (name, d, fs, os);
  return fptr;
}

hpp::constraints::ConvexShapeContactPtr_t createConvexShapeContact_triangles (hpp::pinocchio::DevicePtr_t d, hpp::pinocchio::JointPtr_t j, std::string name)
{
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
  std::vector <vector3_t> square(4);
  square[0] = vector3_t ( 5, 5,0); square[1] = vector3_t ( 5,-5,0);
  square[2] = vector3_t ( 0,-5,0); square[3] = vector3_t ( 0, 5,0);

  std::vector <vector3_t> penta(5);
  penta[0] = vector3_t ( 0, 5,0); penta[1] = vector3_t ( 0,-5,0);
  penta[2] = vector3_t (-2,-6,0); penta[3] = vector3_t (-5, 0,0);
  penta[4] = vector3_t (-2, 6,0);

  std::vector <vector3_t> trapeze(4);
  trapeze[0] = vector3_t (-0.1, 0.1,0); trapeze[1] = vector3_t ( 0.1, 0.1,0);
  trapeze[2] = vector3_t ( 0.2,-0.1,0); trapeze[3] = vector3_t (-0.1,-0.1,0);

  JointAndShapes_t os; os.push_back(JointAndShape_t(j, trapeze));
  JointAndShapes_t fs;
  fs.push_back(JointAndShape_t(JointPtr_t(),square));
  fs.push_back(JointAndShape_t(JointPtr_t(),penta));
  ConvexShapeContactPtr_t fptr = ConvexShapeContact::create (name, d, fs, os);
  return fptr;
}
