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

#define BOOST_TEST_MODULE values
#include <boost/test/included/unit_test.hpp>
#include <boost/assign/list_of.hpp>

#include <hpp/model/device.hh>
#include <hpp/pinocchio/device.hh>

#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/hpp-model/conversions.hh>
#include <hpp/pinocchio/hpp-model/model-loader.hh>

#include "hpp/_constraints/generic-transformation.hh"
#include "hpp/constraints/generic-transformation.hh"
#include "hpp/_constraints/relative-com.hh"
#include "hpp/constraints/relative-com.hh"
#include "hpp/_constraints/com-between-feet.hh"
#include "hpp/constraints/com-between-feet.hh"
#include "hpp/_constraints/symbolic-function.hh"
#include "hpp/constraints/symbolic-function.hh"
#include "hpp/_constraints/distance-between-bodies.hh"
#include "hpp/constraints/distance-between-bodies.hh"
#include "hpp/_constraints/distance-between-points-in-bodies.hh"
#include "hpp/constraints/distance-between-points-in-bodies.hh"
#include <hpp/_constraints/configuration-constraint.hh>
#include <hpp/constraints/configuration-constraint.hh>

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

BOOST_AUTO_TEST_CASE (absolute) {
  model::DevicePtr_t rm;
  pinoc::DevicePtr_t rp;
  setupRobots(rm, rp);

  _c::JointPtr_t eeM = rm->getJointByName ("RWristPitch");
  c ::JointPtr_t eeP = rp->getJointByName ("RWristPitch");

  _c::Transform3f tfM (eeM->currentTransformation ());
  c ::Transform3f tfP (eeP->currentTransformation ());

  c ::Transform3f randP = se3::SE3::Random();
  _c::Transform3f randM (p2m::SE3(randP));

  // This two frames are the position to be compared.
  _c::Transform3f frameM = eeM->linkInJointFrame();
  c ::Transform3f frameP = rp->model().getFramePlacement(eeM->linkName());

  BOOST_REQUIRE(m2p::SE3(tfM * frameM).isApprox(tfP * frameP));

  /*********************** Position **************************/
  // /*
  // Position of the center in world frame.
  check_consistent (rm, rp,
        _c::Position::create ("ModelPosition", rm, eeM, tIdM),
        c ::Position::create ("PinocPosition", rp, eeP, tIdP),
        Compare());

  // Position of the center in some frame.
  check_consistent (rm, rp,
        _c::Position::create ("ModelPosition", rm, eeM, tIdM, randM),
        c ::Position::create ("PinocPosition", rp, eeP, tIdP, randP),
        Compare());

  // Position of a point in some frame.
  check_consistent (rm, rp,
        _c::Position::create ("ModelPosition", rm, eeM, frameM, randM),
        c ::Position::create ("PinocPosition", rp, eeP, frameP, randP),
        Compare());
  // */

  /*********************** Orientation **************************/
  // /*
  // Orientation of a frame in joint frame wrt world frame.
  check_consistent (rm, rp,
        _c::Orientation::create ("ModelOrientation", rm, eeM, frameM, tIdM),
        c ::Orientation::create ("PinocOrientation", rp, eeP, frameP, tIdP),
        Compare());

  // Orientation of a joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::Orientation::create ("ModelOrientation", rm, eeM, frameM * p2m::SE3(frameP.inverse()), randM),
        c ::Orientation::create ("PinocOrientation", rp, eeP, tIdP                               , randP),
        Compare());

  // Orientation of a frame in joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::Orientation::create ("ModelOrientation", rm, eeM, frameM * randM, randM),
        c ::Orientation::create ("PinocOrientation", rp, eeP, frameP * randP, randP),
        Compare());
  // */

  /*********************** Transformation **************************/
  // /*
  // Transformation of a frame in joint frame wrt world frame.
  check_consistent (rm, rp,
        _c::Transformation::create ("ModelTransformation", rm, eeM, frameM, tIdM),
        c ::Transformation::create ("PinocTransformation", rp, eeP, frameP, tIdP),
        Compare());

  // Transformation of a joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::Transformation::create ("ModelTransformation", rm, eeM, frameM * p2m::SE3(frameP.inverse()), randM),
        c ::Transformation::create ("PinocTransformation", rp, eeP, tIdP                               , randP),
        Compare());

  // Transformation of a frame in joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::Transformation::create ("ModelTransformation", rm, eeM, frameM * randM, randM),
        c ::Transformation::create ("PinocTransformation", rp, eeP, frameP * randP, randP),
        Compare());
  // */
}

BOOST_AUTO_TEST_CASE (relative) {
  model::DevicePtr_t rm;
  pinoc::DevicePtr_t rp;
  setupRobots(rm, rp);

  _c::JointPtr_t eeM1 = rm->getJointByName ("RWristPitch");
  c ::JointPtr_t eeP1 = rp->getJointByName ("RWristPitch");

  _c::JointPtr_t eeM2 = rm->getJointByName ("LWristPitch");
  c ::JointPtr_t eeP2 = rp->getJointByName ("LWristPitch");

  _c::Transform3f tfM1 (eeM1->currentTransformation ());
  c ::Transform3f tfP1 (eeP1->currentTransformation ());

  _c::Transform3f tfM2 (eeM2->currentTransformation ());
  c ::Transform3f tfP2 (eeP2->currentTransformation ());

  c ::Transform3f randP1 = se3::SE3::Random();
  _c::Transform3f randM1 (p2m::SE3(randP1));

  c ::Transform3f randP2 = se3::SE3::Random();
  _c::Transform3f randM2 (p2m::SE3(randP2));

  // This two frames are the position to be compared.
  _c::Transform3f frameM1 = eeM1->linkInJointFrame();
  c ::Transform3f frameP1 = rp->model().getFramePlacement(eeM1->linkName());

  _c::Transform3f frameM2 = eeM2->linkInJointFrame();
  c ::Transform3f frameP2 = rp->model().getFramePlacement(eeM2->linkName());

  _c::Transform3f Fp2m1 = frameM1 * p2m::SE3(frameP1.inverse());
  _c::Transform3f Fp2m2 = frameM2 * p2m::SE3(frameP2.inverse());

  BOOST_REQUIRE(m2p::SE3(tfM1 * frameM1).isApprox(tfP1 * frameP1));
  BOOST_REQUIRE(m2p::SE3(tfM2 * frameM2).isApprox(tfP2 * frameP2));

  /*********************** Position **************************/
  // /*
  // Position of a point in some frame.
  check_consistent (rm, rp,
        _c::RelativePosition::create ("ModelRelativePosition", rm, eeM1, eeM2, Fp2m1, Fp2m2),
        c ::RelativePosition::create ("PinocRelativePosition", rp, eeP1, eeP2, tIdP , tIdP ),
        Compare());

  // Position of the center in world frame.
  check_consistent (rm, rp,
        _c::RelativePosition::create ("ModelRelativePosition", rm, eeM1, eeM2, Fp2m1, frameM2 * randM2),
        c ::RelativePosition::create ("PinocRelativePosition", rp, eeP1, eeP2, tIdP , frameP2 * randP2),
        Compare());

  // Position of the center in some frame.
  check_consistent (rm, rp,
        _c::RelativePosition::create ("ModelRelativePosition", rm, eeM1, eeM2, frameM1 * randM1, Fp2m2),
        c ::RelativePosition::create ("PinocRelativePosition", rp, eeP1, eeP2, frameP1 * randP1, tIdP ),
        Compare());

  // Position of a point in some frame.
  check_consistent (rm, rp,
        _c::RelativePosition::create ("ModelRelativePosition", rm, eeM1, eeM2, frameM1 * randM1, frameM2 * randM2),
        c ::RelativePosition::create ("PinocRelativePosition", rp, eeP1, eeP2, frameP1 * randP1, frameP2 * randP2),
        Compare());
  // */

  /*********************** Orientation **************************/
  // /*
  // Orientation of a frame in joint frame wrt world frame.
  check_consistent (rm, rp,
        _c::RelativeOrientation::create ("ModelRelativeOrientation", rm, eeM1, eeM2, Fp2m1, Fp2m2),
        c ::RelativeOrientation::create ("PinocRelativeOrientation", rp, eeP1, eeP2, tIdP , tIdP ),
        Compare());

  // Orientation of a joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::RelativeOrientation::create ("ModelRelativeOrientation", rm, eeM1, eeM2, Fp2m1, frameM2 * randM2),
        c ::RelativeOrientation::create ("PinocRelativeOrientation", rp, eeP1, eeP2, tIdP , frameP2 * randP2),
        Compare());

  // Orientation of a frame in joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::RelativeOrientation::create ("ModelRelativeOrientation", rm, eeM1, eeM2, frameM1 * randM1, Fp2m2),
        c ::RelativeOrientation::create ("PinocRelativeOrientation", rp, eeP1, eeP2, frameP1 * randP1, tIdP ),
        Compare());

  // Orientation of a frame in joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::RelativeOrientation::create ("ModelRelativeOrientation", rm, eeM1, eeM2, frameM1 * randM1, frameM2 * randM2),
        c ::RelativeOrientation::create ("PinocRelativeOrientation", rp, eeP1, eeP2, frameP1 * randP1, frameP2 * randP2),
        Compare());
  // */

  /*********************** Transformation **************************/
  // /*
  // Transformation of a frame in joint frame wrt world frame.
  check_consistent (rm, rp,
        _c::RelativeTransformation::create ("ModelRelativeTransformation", rm, eeM1, eeM2, Fp2m1, Fp2m2),
        c ::RelativeTransformation::create ("PinocRelativeTransformation", rp, eeP1, eeP2, tIdP , tIdP ),
        Compare());

  // Transformation of a joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::RelativeTransformation::create ("ModelRelativeTransformation", rm, eeM1, eeM2, Fp2m1, frameM2 * randM2),
        c ::RelativeTransformation::create ("PinocRelativeTransformation", rp, eeP1, eeP2, tIdP , frameP2 * randP2),
        Compare());

  // Transformation of a frame in joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::RelativeTransformation::create ("ModelRelativeTransformation", rm, eeM1, eeM2, frameM1 * randM1, Fp2m2),
        c ::RelativeTransformation::create ("PinocRelativeTransformation", rp, eeP1, eeP2, frameP1 * randP1, tIdP ),
        Compare());

  // Transformation of a frame in joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::RelativeTransformation::create ("ModelRelativeTransformation", rm, eeM1, eeM2, frameM1 * randM1, frameM2 * randM2),
        c ::RelativeTransformation::create ("PinocRelativeTransformation", rp, eeP1, eeP2, frameP1 * randP1, frameP2 * randP2),
        Compare());
  // */
}

BOOST_AUTO_TEST_CASE (com) {
  model::DevicePtr_t rm;
  pinoc::DevicePtr_t rp;
  setupRobots(rm, rp);

  _c::JointPtr_t eeMR = rm->getJointByName ("RAnkleRoll");
  c ::JointPtr_t eePR = rp->getJointByName ("RAnkleRoll");
  _c::JointPtr_t eeML = rm->getJointByName ("LAnkleRoll");
  c ::JointPtr_t eePL = rp->getJointByName ("LAnkleRoll");

  _c::Transform3f tfMR (eeMR->currentTransformation ());
  c ::Transform3f tfPR (eePR->currentTransformation ());

  _c::vector3_t targetMR = tfMR.getRotation() * (rm->positionCenterOfMass() - tfMR.getTranslation());
  c ::vector3_t targetPR = tfPR.actInv(rm->positionCenterOfMass().derived());

  // This two frames are the position to be compared.
  _c::Transform3f frameMR = eeMR->linkInJointFrame();
  c ::Transform3f framePR = rp->model().getFramePlacement(eeMR->linkName());

  BOOST_REQUIRE(m2p::SE3(tfMR * frameMR).isApprox(tfPR * framePR));

  _c::CenterOfMassComputationPtr_t comM = _c::CenterOfMassComputation::create(rm);
  c ::CenterOfMassComputationPtr_t comP = c ::CenterOfMassComputation::create(rp);
  comM->add(rm->rootJoint());
  comP->add(rp->rootJoint());
  // FIXME This does not work for an unknown reason
  // comM->add(rm->getJointByName("LHipYaw"));
  // comP->add(rp->getJointByName("LHipYaw"));
  comM->computeMass();
  comP->computeMass();

  comM->compute(_c::Device::COM);
  comP->compute(c ::Device::COM);
  BOOST_CHECK(comM->com().isApprox(comP->com()));
  comM->compute(_c::Device::ALL);
  comP->compute(c ::Device::ALL);
  BOOST_CHECK(comM->com().isApprox(comP->com()));
  BOOST_CHECK(comM->jacobian().isApprox(
        comP->jacobian() * m2p::Xq(rp->rootJoint()->currentTransformation())));

  /*********************** COM via SymbolicFunction **************************/
  // /*
  typedef _c::SymbolicFunction<_c::PointCom> MCom;
  typedef c ::SymbolicFunction<c ::PointCom> PCom;

  check_consistent (rm, rp,
        MCom::create ("ModelPointCom", rm, _c::PointCom::create(comM)),
        PCom::create ("PinocPointCom", rp, c ::PointCom::create(comP)),
        Compare());
  // */

  /*********************** Relative COM **************************/
  // /*
  check_consistent (rm, rp,
        _c::RelativeCom::create (rm, comM, eeMR, targetMR),
        c ::RelativeCom::create (rp, comP, eePR, targetPR),
        Compare());
  // */

  /*********************** COM between feet **************************/
  // /*
  check_consistent (rm, rp,
        _c::ComBetweenFeet::create ("ModelComBetweenFeet", rm, eeML, eeMR, tIdM.getTranslation(), tIdM.getTranslation(), eeMR, tIdM.getTranslation()),
        c ::ComBetweenFeet::create ("PinocComBetweenFeet", rp, eePL, eePR, tIdP.translation()   , tIdP.translation()   , eePR, tIdP.translation()),
        Compare());
  // */
}

BOOST_AUTO_TEST_CASE (distance) {
  model::DevicePtr_t rm;
  pinoc::DevicePtr_t rp;
  setupRobots(rm, rp, true);

  _c::JointPtr_t eeMR = rm->getJointByName ("RAnkleRoll");
  c ::JointPtr_t eePR = rp->getJointByName ("RAnkleRoll");
  _c::JointPtr_t eeML = rm->getJointByName ("LAnkleRoll");
  c ::JointPtr_t eePL = rp->getJointByName ("LAnkleRoll");

  c ::Transform3f randPR = se3::SE3::Random();
  _c::Transform3f randMR (p2m::SE3(randPR));

  c ::Transform3f randPL = se3::SE3::Random();
  _c::Transform3f randML (p2m::SE3(randPL));

  // This two frames are the position to be compared.
  _c::Transform3f frameMR = eeMR->linkInJointFrame();
  c ::Transform3f framePR = rp->model().getFramePlacement(eeMR->linkName());

  _c::Transform3f frameML = eeML->linkInJointFrame();
  c ::Transform3f framePL = rp->model().getFramePlacement(eeML->linkName());

  /*********************** Distance between bodies **************************/
  // /*
  check_consistent (rm, rp,
        _c::DistanceBetweenBodies::create ("ModelDistanceBetweenBodies", rm, eeMR, eeML),
        c ::DistanceBetweenBodies::create ("PinocDistanceBetweenBodies", rp, eePR, eePL),
        Compare());
  // */

  /*********************** Distance between point in bodies **************************/
  // /*
  check_consistent (rm, rp,
        _c::DistanceBetweenPointsInBodies::create ("ModelDistanceBetweenPointInBodies", rm, eeMR, eeML, (frameMR * randMR).getTranslation(), (frameML * randML).getTranslation()),
        c ::DistanceBetweenPointsInBodies::create ("PinocDistanceBetweenPointInBodies", rp, eePR, eePL, (framePR * randPR).translation(), (framePL * randPL).translation()),
        Compare());
  // */
}

// ConfigurationConstraint
BOOST_AUTO_TEST_CASE (others) {
  model::DevicePtr_t rm;
  pinoc::DevicePtr_t rp;
  setupRobots(rm, rp, true);

  _c::Configuration_t goalM = rm->neutralConfiguration();
  c ::Configuration_t goalP = m2p::q(goalM);

  /*********************** ConfigurationConstraint **************************/
  // /*
  check_consistent (rm, rp,
        _c::ConfigurationConstraint::create ("Model ConfigurationConstraint", rm, goalM),
        c ::ConfigurationConstraint::create ("Pinoc ConfigurationConstraint", rp, goalP),
        Compare());
  // */
}
