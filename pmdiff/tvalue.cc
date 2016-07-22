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
#include "hpp/_constraints/relative-com.hh"
#include "hpp/constraints/generic-transformation.hh"
#include "hpp/constraints/relative-com.hh"


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

struct TransformationCompare {
  // const matrix3_t R;
  // const size_type rowTr, sizeTr;
  // TransformationCompare (const matrix3_t _R, size_type _rowTr, size_type _sizeTr)
    // : R (_R), rowTr(_rowTr), sizeTr(_sizeTr) {}

  // void value(const vector_t& value1, const vector_t& value2) const {
    // vector_t dtr  = R.middleRows(rowTr, sizeTr) * value2.segment(rowTr,sizeTr) + value1.segment(rowTr,sizeTr);
    // BOOST_CHECK_MESSAGE(dtr.isZero(1e-10),
			// "Value not matching. Norm of value1 is "
			// << value1.segment(rowTr,sizeTr).norm () <<  ", norm of value2 is "
			// << value2.segment(rowTr,sizeTr).norm () << ", norm of error is "
			// << dtr.norm());
    // //matrix3_t A = exponential(value1.tail<3>()) * transpose(R);
    // //matrix3_t B = exponential(value2.tail<3>()) * R;
    // matrix3_t A = exponential(value1.tail<3>());
    // matrix3_t B = exponential(value2.tail<3>());
    // matrix_t diff = A.transpose() - B;
    // // std::cout << diff << std::endl;
    // BOOST_CHECK_MESSAGE(diff.isZero(1e-10),
			// "Value not matching. Norm of error is "
			// << diff.norm());
  // }

  // void jacobian(const matrix_t& jacobian1, const matrix_t& jacobian2) const {
    // matrix_t diffJtr = R.middleRows(rowTr, sizeTr) * jacobian2.middleRows(rowTr, sizeTr) + jacobian1.middleRows(rowTr, sizeTr);
    // //matrix_t diffJrot = jacobian2 - jacobian1;
    // // std::cout << diffJ.norm() << std::endl;
    // BOOST_CHECK_MESSAGE(diffJtr.isZero(1e-10),
			// "Jacobian not matching. Norm of J1 is "
			// << jacobian1.middleRows(rowTr, sizeTr).norm () << ", norm of J2 is "
			// << jacobian2.middleRows(rowTr, sizeTr).norm () << ", norm of error is "
			// << diffJtr.norm());
    // std::cout << "Rotation part of jacobian not checked." << std::endl;
  // }
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

BOOST_AUTO_TEST_SUITE ( PinocchioModelValue )

BOOST_AUTO_TEST_CASE (absolute) {
  model::DevicePtr_t rm = hppModel();
  pinoc::DevicePtr_t rp = hppPinocchio();

  _c::Configuration_t qm = rm->neutralConfiguration();
  // _c::Configuration_t q = rp->neutralConfiguration();
  c ::Configuration_t qp = m2p::q(qm);
  rm->currentConfiguration(qm); rm->computeForwardKinematics();
  rp->currentConfiguration(qp); rp->computeForwardKinematics();

  /// Set root joint bound.
  rm->rootJoint()->lowerBound(0,-1); rm->rootJoint()->lowerBound(1,-1); rm->rootJoint()->lowerBound(2,-1);
  rm->rootJoint()->upperBound(0, 1); rm->rootJoint()->upperBound(1, 1); rm->rootJoint()->upperBound(2, 1);
  rp->rootJoint()->lowerBound(0,-1); rp->rootJoint()->lowerBound(1,-1); rp->rootJoint()->lowerBound(2,-1);
  rp->rootJoint()->upperBound(0, 1); rp->rootJoint()->upperBound(1, 1); rp->rootJoint()->upperBound(2, 1);

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
  /*
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
  /*
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
  /*
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
  model::DevicePtr_t rm = hppModel();
  pinoc::DevicePtr_t rp = hppPinocchio();

  _c::Configuration_t qm = rm->neutralConfiguration();
  // _c::Configuration_t q = rp->neutralConfiguration();
  c ::Configuration_t qp = m2p::q(qm);
  rm->currentConfiguration(qm); rm->computeForwardKinematics();
  rp->currentConfiguration(qp); rp->computeForwardKinematics();

  /// Set root joint bound.
  rm->rootJoint()->lowerBound(0,-1); rm->rootJoint()->lowerBound(1,-1); rm->rootJoint()->lowerBound(2,-1);
  rm->rootJoint()->upperBound(0, 1); rm->rootJoint()->upperBound(1, 1); rm->rootJoint()->upperBound(2, 1);
  rp->rootJoint()->lowerBound(0,-1); rp->rootJoint()->lowerBound(1,-1); rp->rootJoint()->lowerBound(2,-1);
  rp->rootJoint()->upperBound(0, 1); rp->rootJoint()->upperBound(1, 1); rp->rootJoint()->upperBound(2, 1);

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
  model::DevicePtr_t rm = hppModel();
  pinoc::DevicePtr_t rp = hppPinocchio();

  _c::Configuration_t qm = rm->neutralConfiguration();
  // _c::Configuration_t q = rp->neutralConfiguration();
  c ::Configuration_t qp = m2p::q(qm);
  rm->controlComputation((_c::Device::Computation_t)(_c::Device::COM | _c::Device::JACOBIAN | _c::Device::JOINT_POSITION));
  rp->controlComputation(( c::Device::Computation_t)( c::Device::COM |  c::Device::JACOBIAN |  c::Device::JOINT_POSITION));
  rm->currentConfiguration(qm); rm->computeForwardKinematics();
  rp->currentConfiguration(qp); rp->computeForwardKinematics();

  /// Set root joint bound.
  rm->rootJoint()->lowerBound(0,-1); rm->rootJoint()->lowerBound(1,-1); rm->rootJoint()->lowerBound(2,-1);
  rm->rootJoint()->upperBound(0, 1); rm->rootJoint()->upperBound(1, 1); rm->rootJoint()->upperBound(2, 1);
  rp->rootJoint()->lowerBound(0,-1); rp->rootJoint()->lowerBound(1,-1); rp->rootJoint()->lowerBound(2,-1);
  rp->rootJoint()->upperBound(0, 1); rp->rootJoint()->upperBound(1, 1); rp->rootJoint()->upperBound(2, 1);

  _c::JointPtr_t eeM = rm->getJointByName ("RAnkleRoll");
  c ::JointPtr_t eeP = rp->getJointByName ("RAnkleRoll");

  _c::Transform3f tfM (eeM->currentTransformation ());
  c ::Transform3f tfP (eeP->currentTransformation ());

  _c::vector3_t targetM = tfM.getRotation() * (rm->positionCenterOfMass() - tfM.getTranslation());
  c ::vector3_t targetP = tfP.actInv(rm->positionCenterOfMass().derived());

  // This two frames are the position to be compared.
  _c::Transform3f frameM = eeM->linkInJointFrame();
  c ::Transform3f frameP = rp->model()->getFramePlacement(eeM->linkName());

  BOOST_REQUIRE(m2p::SE3(tfM * frameM).isApprox(tfP * frameP));

  /*********************** Relative COM **************************/
  // /*
  check_consistent (rm, rp,
        _c::RelativeCom::create (rm, eeM, targetM),
        c ::RelativeCom::create (rp, eeP, targetP),
        ProportionalCompare(1));
  // */
}

BOOST_AUTO_TEST_SUITE_END ()
