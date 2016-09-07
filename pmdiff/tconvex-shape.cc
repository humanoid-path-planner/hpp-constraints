// Copyright (c) 2016, Joseph Mirabel
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

#define BOOST_TEST_MODULE convex-shape
#include <boost/test/included/unit_test.hpp>
#include <boost/assign/list_of.hpp>

#include <hpp/model/device.hh>
#include <hpp/pinocchio/device.hh>

#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/hpp-model/conversions.hh>
#include <hpp/pinocchio/hpp-model/model-loader.hh>

#include "hpp/_constraints/convex-shape-contact.hh"
#include "hpp/constraints/convex-shape-contact.hh"

#include <stdlib.h>
#include <limits>
#include <math.h>

const static size_t NUMBER_RANDOM_SAMPLES = 10;
const bool verbose = false;
const bool verboseNum = false;

using std::numeric_limits;
using boost::assign::list_of;

typedef std::vector<bool> BoolVector_t;

namespace c  = hpp::constraints;
namespace _c = hpp::_constraints;

namespace model = hpp::model    ;
namespace pinoc = hpp::pinocchio;

#include <../pmdiff/tools.cc>

_c::Transform3f tIdM = _c::Transform3f();
c ::Transform3f tIdP = c ::Transform3f::Identity();

struct Shapes {
  typedef std::vector<fcl::Vec3f> FclPoints_t;
  typedef std::vector<FclPoints_t> FclPointss_t;
  typedef std::vector<c::vector3_t> Points_t;
  typedef std::vector<Points_t> Pointss_t;
  FclPointss_t fclObject;
  FclPointss_t fclFloor;
/*
  template <typename In, typename Out> void set(const In& in, Out& out) { out = in; }

  template <typename In, typename Out>
  std::vector<Out> set(const std::vector<In>& in, std::vector<Out>& out) {
    out.resize(in.size());
    for (std::size_t i = 0; i < in.size(); ++i) set<In,Out>(in[i], out[i]);
  }

  template <typename In, typename Out>
  std::vector<Out> set(const std::vector<In>& in, std::vector<Out>& out) {
    out.resize(in.size());
    for (std::size_t i = 0; i < in.size(); ++i) set<In,Out>(in[i], out[i]);
  }

  Pointss_t object() {
    Pointss_t pts;
    set(fclObject, pts);
    return pts;
  }

  Pointss_t floor() {
    Pointss_t pts;
    set(fclObject, pts);
    return pts;
  }
//*/
  template <bool FLOOR> void add (_c::ConvexShapeContactPtr_t csc, _c::JointPtr_t j)
  {
    FclPointss_t& ptss = (FLOOR ? fclFloor : fclObject);
    for (std::size_t i = 0; i < ptss.size (); ++i)
    {
      if (FLOOR) csc->addFloor (_c::ConvexShape (ptss[i], j));
      else       csc->addObject(_c::ConvexShape (ptss[i], j));
    }
  }

  template <bool FLOOR> void add (c ::ConvexShapeContactPtr_t csc, c ::JointPtr_t j, c ::Transform3f M)
  {
    FclPointss_t& ptss = (FLOOR ? fclFloor : fclObject);
    for (std::size_t i = 0; i < ptss.size (); ++i)
    {
      Points_t pts(ptss[i].size());
      for (std::size_t k = 0; k < ptss[i].size(); ++k) pts[k] = M.act(ptss[i][k].derived());
      if (FLOOR) csc->addFloor (c::ConvexShape (pts, j));
      else       csc->addObject(c::ConvexShape (pts, j));
    }
  }
};

/***************** Punctual ************************/
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
Shapes convexShapes_punctual ()
{
  Shapes shapes;
  std::vector <fcl::Vec3f> square(4);
  square[0] = fcl::Vec3f ( 5, 5,0); square[1] = fcl::Vec3f ( 5,-5,0);
  square[2] = fcl::Vec3f ( 0,-5,0); square[3] = fcl::Vec3f ( 0, 5,0);
  shapes.fclFloor.push_back(square);

  std::vector <fcl::Vec3f> penta(5);
  penta[0] = fcl::Vec3f ( 0, 5,0); penta[1] = fcl::Vec3f ( 0,-5,0);
  penta[2] = fcl::Vec3f (-2,-6,0); penta[3] = fcl::Vec3f (-5, 0,0);
  penta[4] = fcl::Vec3f (-2, 6,0);
  shapes.fclFloor.push_back(penta);

  std::vector <fcl::Vec3f> point(1, fcl::Vec3f (0,0,0));
  shapes.fclObject.push_back(point);

  return shapes;
}

/***************** Convex ************************/
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
Shapes convexShapes_convex ()
{
  Shapes shapes;
  std::vector <fcl::Vec3f> square(4);
  square[0] = fcl::Vec3f ( 5, 5,0); square[1] = fcl::Vec3f ( 5,-5,0);
  square[2] = fcl::Vec3f ( 0,-5,0); square[3] = fcl::Vec3f ( 0, 5,0);
  shapes.fclFloor.push_back(square);

  std::vector <fcl::Vec3f> penta(5);
  penta[0] = fcl::Vec3f ( 0, 5,0); penta[1] = fcl::Vec3f ( 0,-5,0);
  penta[2] = fcl::Vec3f (-2,-6,0); penta[3] = fcl::Vec3f (-5, 0,0);
  penta[4] = fcl::Vec3f (-2, 6,0);
  shapes.fclFloor.push_back(penta);

  std::vector <fcl::Vec3f> trapeze(4);
  trapeze[0] = fcl::Vec3f (-0.1, 0.1,0); trapeze[1] = fcl::Vec3f ( 0.1, 0.1,0);
  trapeze[2] = fcl::Vec3f ( 0.2,-0.1,0); trapeze[3] = fcl::Vec3f (-0.1,-0.1,0);
  shapes.fclObject.push_back(trapeze);

  return shapes;
}

/***************** Triangle ************************/
Shapes convexShapes_triangle ()
{
  fcl::Vec3f x (1,0,0), y (0,1,0), z (0,0,1);
  fcl::Vec3f p[12];
  p[0] = fcl::Vec3f (-5,-5,0); p[1] = fcl::Vec3f (-5, 5,0); p[2] = fcl::Vec3f ( 5,-5,0);
  p[3] = fcl::Vec3f ( 5, 5,0); p[4] = fcl::Vec3f (-5, 5,0); p[5] = fcl::Vec3f ( 5,-5,0);
  p[6] = fcl::Vec3f ( 0, 0,1); p[7] = fcl::Vec3f (  1,0,1); p[8] = fcl::Vec3f (0,  1,1);
  p[9] = fcl::Vec3f ( 0, 0,0); p[10] = fcl::Vec3f (0.1,0,0); p[11] = fcl::Vec3f (0,0.1,0);
  Shapes shapes;
  shapes.fclObject.push_back(list_of(p[0])(p[1])(p[2]));
  shapes.fclFloor.push_back(list_of(p[3])(p[4])(p[5]));
  shapes.fclFloor.push_back(list_of(p[6])(p[7])(p[8]));
  shapes.fclFloor.push_back(list_of(p[9])(p[10])(p[11]));
  return shapes;
}

BOOST_AUTO_TEST_CASE (convexShapeContact) {
  model::DevicePtr_t rm;
  pinoc::DevicePtr_t rp;
  setupRobots(rm, rp);

  _c::JointPtr_t eeM = rm->getJointByName ("RWristPitch");
  c ::JointPtr_t eeP = rp->getJointByName ("RWristPitch");

  BOOST_CHECK (rp->model().existFrame(eeM->linkName(), se3::BODY));
  _c::Transform3f frameM = eeM->linkInJointFrame();
  c ::Transform3f frameP = rp->model().frames[rp->model().getFrameId(eeM->linkName(), se3::BODY)].placement;

  c ::Transform3f Fm2p = frameP * m2p::SE3(frameM).inverse();

  _c::ConvexShapeContactPtr_t cscM;
  c ::ConvexShapeContactPtr_t cscP;

  /*********************** Convex **************************/
  // /*
  cscM = _c::ConvexShapeContact::create(rm);
  cscP = c ::ConvexShapeContact::create(rp);

  Shapes convex = convexShapes_convex ();
  convex.add<true> (cscM, _c::JointPtr_t());
  convex.add<false>(cscM, eeM);
  convex.add<true> (cscP,  c::JointPtr_t(), tIdP);
  convex.add<false>(cscP, eeP, Fm2p);

  check_consistent (rm, rp, cscM, cscP, Compare());
  // */

  /*********************** Triangle **************************/
  // /*
  cscM = _c::ConvexShapeContact::create(rm);
  cscP = c ::ConvexShapeContact::create(rp);

  Shapes tri = convexShapes_triangle ();
  tri.add<true> (cscM, _c::JointPtr_t());
  tri.add<false>(cscM, eeM);
  tri.add<true> (cscP,  c::JointPtr_t(), tIdP);
  tri.add<false>(cscP, eeP, Fm2p);

  check_consistent (rm, rp, cscM, cscP, Compare());
  // */

  /*********************** Triangle **************************/
  // /*
  cscM = _c::ConvexShapeContact::create(rm);
  cscP = c ::ConvexShapeContact::create(rp);

  Shapes punctual = convexShapes_punctual ();
  punctual.add<true> (cscM, _c::JointPtr_t());
  punctual.add<false>(cscM, eeM);
  punctual.add<true> (cscP,  c::JointPtr_t(), tIdP);
  punctual.add<false>(cscP, eeP, Fm2p);

  check_consistent (rm, rp, cscM, cscP, Compare());
  // */
}
