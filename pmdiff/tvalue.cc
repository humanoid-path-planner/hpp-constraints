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

_c::Transform3f tIdM = _c::Transform3f();
c ::Transform3f tIdP = c ::Transform3f::Identity();

void setupRobots(model::DevicePtr_t& rm, pinoc::DevicePtr_t& rp, bool geom = false)
{
  rm = hppModel();
  rp = hppPinocchio(geom);

  _c::Configuration_t qm = rm->neutralConfiguration();
  // _c::Configuration_t q = rp->neutralConfiguration();
  c ::Configuration_t qp = m2p::q(qm);
  if (geom) {
    rm->controlComputation((_c::Device::Computation_t)(                       _c::Device::COM | _c::Device::JACOBIAN | _c::Device::JOINT_POSITION));
    rp->controlComputation(( c::Device::Computation_t)( c::Device::GEOMETRY |  c::Device::COM |  c::Device::JACOBIAN |  c::Device::JOINT_POSITION));
  } else {
    rm->controlComputation((_c::Device::Computation_t)(_c::Device::COM | _c::Device::JACOBIAN | _c::Device::JOINT_POSITION));
    rp->controlComputation(( c::Device::Computation_t)( c::Device::COM |  c::Device::JACOBIAN |  c::Device::JOINT_POSITION));
  }
  rm->currentConfiguration(qm); rm->computeForwardKinematics();
  rp->currentConfiguration(qp); rp->computeForwardKinematics();

  /// Set root joint bound.
  rm->rootJoint()->lowerBound(0,-1); rm->rootJoint()->lowerBound(1,-1); rm->rootJoint()->lowerBound(2,-1);
  rm->rootJoint()->upperBound(0, 1); rm->rootJoint()->upperBound(1, 1); rm->rootJoint()->upperBound(2, 1);
  rp->rootJoint()->lowerBound(0,-1); rp->rootJoint()->lowerBound(1,-1); rp->rootJoint()->lowerBound(2,-1);
  rp->rootJoint()->upperBound(0, 1); rp->rootJoint()->upperBound(1, 1); rp->rootJoint()->upperBound(2, 1);
}

struct ProportionalCompare {
  const c::value_type alpha;
  ProportionalCompare (c::value_type _alpha = 1) : alpha (_alpha) {}

  void value(const _c::vector_t& valueM, const c::vector_t& valueP) const {
    c::vector_t d = valueM - alpha * valueP;
    if (verboseNum) {
      std::cout << "---- hpp::model ------" << std::endl;
      std::cout << valueM.transpose() << std::endl;
      std::cout << "---- hpp::pinocchio ------" << std::endl;
      std::cout << valueP.transpose() << std::endl;
      std::cout << "---- model - " << alpha << " * piniocchio ------" << std::endl;
      std::cout << d.transpose() << std::endl;
    }
    BOOST_CHECK_MESSAGE(d.isZero(1e-10),
			"Value not matching. Norm of value from model is "
			<< valueM.norm () <<  ", norm of value Pinocchio is "
			<< valueP.norm () << ", norm of error is "
			<< d.norm());
  }

  void jacobian(const _c::matrix_t& jacobianM, const c::matrix_t& jacobianP) const {
    c::matrix_t diffJ = jacobianM - jacobianP;
    if (verboseNum) {
      std::cout << "---- hpp::model ------" << std::endl;
      std::cout << jacobianM.transpose() << std::endl;
      std::cout << "---- hpp::pinocchio ------" << std::endl;
      std::cout << jacobianP.transpose() << std::endl;
      std::cout << "---- model - " << alpha << " * piniocchio ------" << std::endl;
      std::cout << diffJ.transpose() << std::endl;
    }
    // // std::cout << diffJ.norm() << std::endl;
    BOOST_CHECK_MESSAGE(diffJ.isZero(1e-10),
			"Jacobian not matching. Norm of Jacobian from model is "
			<< jacobianM.norm () << ", norm of Jacobian from Pinocchio is "
			<< jacobianP.norm () << ", norm of error is "
			<< diffJ.norm());
  }
};

template <typename Compare>
void check_consistent (_c::DevicePtr_t rm, c::DevicePtr_t rp,
    _c::DifferentiableFunctionPtr_t fm, c::DifferentiableFunctionPtr_t fp,
    const Compare comp = Compare())
{
  // std::cout << f->name() << '\n' << g->name() << '\n';
  BOOST_CHECK(fm->outputSize()==fp->outputSize());
  BOOST_CHECK(fm->inputSize()==fp->inputSize());
  BOOST_CHECK(fm->outputDerivativeSize()==fp->outputDerivativeSize());
  BOOST_CHECK(fm->inputDerivativeSize()==fp->inputDerivativeSize());

  _c::vector_t valueM (fm->outputSize ());
  c ::vector_t valueP (fp->outputSize ());
  _c::matrix_t jacobianM (fm->outputSize (), rm->numberDof ());
  c ::matrix_t jacobianP (fp->outputSize (), rp->numberDof ());
  for (size_t i = 0; i < NUMBER_RANDOM_SAMPLES; i++) {
    c ::Configuration_t qp = se3::randomConfiguration(*rp->model());
    _c::Configuration_t qm = p2m::q(qp);
    (*fm) (valueM, qm);
    (*fp) (valueP, qp);
    comp.value (valueM, valueP);
    fm->jacobian (jacobianM, qm);
    fp->jacobian (jacobianP, qp);
    comp.jacobian (jacobianM, jacobianP * m2p::Xq(rp->rootJoint()->currentTransformation()));
  }
}

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
  c ::Transform3f frameP = rp->model()->getFramePlacement(eeM->linkName());

  BOOST_REQUIRE(m2p::SE3(tfM * frameM).isApprox(tfP * frameP));

  /*********************** Position **************************/
  // /*
  // Position of the center in world frame.
  check_consistent (rm, rp,
        _c::Position::create ("ModelPosition", rm, eeM, tIdM),
        c ::Position::create ("PinocPosition", rp, eeP, tIdP),
        ProportionalCompare(1));

  // Position of the center in some frame.
  check_consistent (rm, rp,
        _c::Position::create ("ModelPosition", rm, eeM, tIdM, randM),
        c ::Position::create ("PinocPosition", rp, eeP, tIdP, randP),
        ProportionalCompare(1));

  // Position of a point in some frame.
  check_consistent (rm, rp,
        _c::Position::create ("ModelPosition", rm, eeM, frameM, randM),
        c ::Position::create ("PinocPosition", rp, eeP, frameP, randP),
        ProportionalCompare(1));
  // */

  /*********************** Orientation **************************/
  // /*
  // Orientation of a frame in joint frame wrt world frame.
  check_consistent (rm, rp,
        _c::Orientation::create ("ModelOrientation", rm, eeM, frameM, tIdM),
        c ::Orientation::create ("PinocOrientation", rp, eeP, frameP, tIdP),
        ProportionalCompare(1));

  // Orientation of a joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::Orientation::create ("ModelOrientation", rm, eeM, frameM * p2m::SE3(frameP.inverse()), randM),
        c ::Orientation::create ("PinocOrientation", rp, eeP, tIdP                               , randP),
        ProportionalCompare(1));

  // Orientation of a frame in joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::Orientation::create ("ModelOrientation", rm, eeM, frameM * randM, randM),
        c ::Orientation::create ("PinocOrientation", rp, eeP, frameP * randP, randP),
        ProportionalCompare(1));
  // */

  /*********************** Transformation **************************/
  // /*
  // Transformation of a frame in joint frame wrt world frame.
  check_consistent (rm, rp,
        _c::Transformation::create ("ModelTransformation", rm, eeM, frameM, tIdM),
        c ::Transformation::create ("PinocTransformation", rp, eeP, frameP, tIdP),
        ProportionalCompare(1));

  // Transformation of a joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::Transformation::create ("ModelTransformation", rm, eeM, frameM * p2m::SE3(frameP.inverse()), randM),
        c ::Transformation::create ("PinocTransformation", rp, eeP, tIdP                               , randP),
        ProportionalCompare(1));

  // Transformation of a frame in joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::Transformation::create ("ModelTransformation", rm, eeM, frameM * randM, randM),
        c ::Transformation::create ("PinocTransformation", rp, eeP, frameP * randP, randP),
        ProportionalCompare(1));
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
  c ::Transform3f frameP1 = rp->model()->getFramePlacement(eeM1->linkName());

  _c::Transform3f frameM2 = eeM2->linkInJointFrame();
  c ::Transform3f frameP2 = rp->model()->getFramePlacement(eeM2->linkName());

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
        ProportionalCompare(1));

  // Position of the center in world frame.
  check_consistent (rm, rp,
        _c::RelativePosition::create ("ModelRelativePosition", rm, eeM1, eeM2, Fp2m1, frameM2 * randM2),
        c ::RelativePosition::create ("PinocRelativePosition", rp, eeP1, eeP2, tIdP , frameP2 * randP2),
        ProportionalCompare(1));

  // Position of the center in some frame.
  check_consistent (rm, rp,
        _c::RelativePosition::create ("ModelRelativePosition", rm, eeM1, eeM2, frameM1 * randM1, Fp2m2),
        c ::RelativePosition::create ("PinocRelativePosition", rp, eeP1, eeP2, frameP1 * randP1, tIdP ),
        ProportionalCompare(1));

  // Position of a point in some frame.
  check_consistent (rm, rp,
        _c::RelativePosition::create ("ModelRelativePosition", rm, eeM1, eeM2, frameM1 * randM1, frameM2 * randM2),
        c ::RelativePosition::create ("PinocRelativePosition", rp, eeP1, eeP2, frameP1 * randP1, frameP2 * randP2),
        ProportionalCompare(1));
  // */

  /*********************** Orientation **************************/
  // /*
  // Orientation of a frame in joint frame wrt world frame.
  check_consistent (rm, rp,
        _c::RelativeOrientation::create ("ModelRelativeOrientation", rm, eeM1, eeM2, Fp2m1, Fp2m2),
        c ::RelativeOrientation::create ("PinocRelativeOrientation", rp, eeP1, eeP2, tIdP , tIdP ),
        ProportionalCompare(1));

  // Orientation of a joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::RelativeOrientation::create ("ModelRelativeOrientation", rm, eeM1, eeM2, Fp2m1, frameM2 * randM2),
        c ::RelativeOrientation::create ("PinocRelativeOrientation", rp, eeP1, eeP2, tIdP , frameP2 * randP2),
        ProportionalCompare(1));

  // Orientation of a frame in joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::RelativeOrientation::create ("ModelRelativeOrientation", rm, eeM1, eeM2, frameM1 * randM1, Fp2m2),
        c ::RelativeOrientation::create ("PinocRelativeOrientation", rp, eeP1, eeP2, frameP1 * randP1, tIdP ),
        ProportionalCompare(1));

  // Orientation of a frame in joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::RelativeOrientation::create ("ModelRelativeOrientation", rm, eeM1, eeM2, frameM1 * randM1, frameM2 * randM2),
        c ::RelativeOrientation::create ("PinocRelativeOrientation", rp, eeP1, eeP2, frameP1 * randP1, frameP2 * randP2),
        ProportionalCompare(1));
  // */

  /*********************** Transformation **************************/
  // /*
  // Transformation of a frame in joint frame wrt world frame.
  check_consistent (rm, rp,
        _c::RelativeTransformation::create ("ModelRelativeTransformation", rm, eeM1, eeM2, Fp2m1, Fp2m2),
        c ::RelativeTransformation::create ("PinocRelativeTransformation", rp, eeP1, eeP2, tIdP , tIdP ),
        ProportionalCompare(1));

  // Transformation of a joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::RelativeTransformation::create ("ModelRelativeTransformation", rm, eeM1, eeM2, Fp2m1, frameM2 * randM2),
        c ::RelativeTransformation::create ("PinocRelativeTransformation", rp, eeP1, eeP2, tIdP , frameP2 * randP2),
        ProportionalCompare(1));

  // Transformation of a frame in joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::RelativeTransformation::create ("ModelRelativeTransformation", rm, eeM1, eeM2, frameM1 * randM1, Fp2m2),
        c ::RelativeTransformation::create ("PinocRelativeTransformation", rp, eeP1, eeP2, frameP1 * randP1, tIdP ),
        ProportionalCompare(1));

  // Transformation of a frame in joint frame wrt to a frame in world frame.
  check_consistent (rm, rp,
        _c::RelativeTransformation::create ("ModelRelativeTransformation", rm, eeM1, eeM2, frameM1 * randM1, frameM2 * randM2),
        c ::RelativeTransformation::create ("PinocRelativeTransformation", rp, eeP1, eeP2, frameP1 * randP1, frameP2 * randP2),
        ProportionalCompare(1));
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
  c ::Transform3f framePR = rp->model()->getFramePlacement(eeMR->linkName());

  BOOST_REQUIRE(m2p::SE3(tfMR * frameMR).isApprox(tfPR * framePR));

  _c::CenterOfMassComputationPtr_t comM = _c::CenterOfMassComputation::create(rm);
  c ::CenterOfMassComputationPtr_t comP = c ::CenterOfMassComputation::create(rp);
  comM->add(rm->rootJoint());
  comP->add(rp->rootJoint());
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
        ProportionalCompare(1));
  // */

  /*********************** Relative COM **************************/
  // /*
  check_consistent (rm, rp,
        _c::RelativeCom::create (rm, comM, eeMR, targetMR),
        c ::RelativeCom::create (rp, comP, eePR, targetPR),
        ProportionalCompare(1));
  // */

  /*********************** COM between feet **************************/
  // /*
  check_consistent (rm, rp,
        _c::ComBetweenFeet::create ("ModelComBetweenFeet", rm, eeML, eeMR, tIdM.getTranslation(), tIdM.getTranslation(), eeMR, tIdM.getTranslation()),
        c ::ComBetweenFeet::create ("PinocComBetweenFeet", rp, eePL, eePR, tIdP.translation()   , tIdP.translation()   , eePR, tIdP.translation()),
        ProportionalCompare(1));
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
  c ::Transform3f framePR = rp->model()->getFramePlacement(eeMR->linkName());

  _c::Transform3f frameML = eeML->linkInJointFrame();
  c ::Transform3f framePL = rp->model()->getFramePlacement(eeML->linkName());

  /*********************** Distance between bodies **************************/
  // /*
  check_consistent (rm, rp,
        _c::DistanceBetweenBodies::create ("ModelDistanceBetweenBodies", rm, eeMR, eeML),
        c ::DistanceBetweenBodies::create ("ModelDistanceBetweenBodies", rp, eePR, eePL),
        ProportionalCompare(1));
  // */

  /*********************** Distance between point in bodies **************************/
  // /*
  check_consistent (rm, rp,
        _c::DistanceBetweenPointsInBodies::create ("ModelDistanceBetweenPointInBodies", rm, eeMR, eeML, (frameMR * randMR).getTranslation(), (frameML * randML).getTranslation()),
        c ::DistanceBetweenPointsInBodies::create ("ModelDistanceBetweenPointInBodies", rp, eePR, eePL, (framePR * randPR).translation(), (framePL * randPL).translation()),
        ProportionalCompare(1));
  // */
}
