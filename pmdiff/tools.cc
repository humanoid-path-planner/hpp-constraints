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

struct Compare {
  void value(const _c::vector_t& valueM, const c::vector_t& valueP) const {
    c::vector_t d = valueM - valueP;
    if (verboseNum) {
      std::cout << "---- hpp::model ------" << std::endl;
      std::cout << valueM.transpose() << std::endl;
      std::cout << "---- hpp::pinocchio ------" << std::endl;
      std::cout << valueP.transpose() << std::endl;
      std::cout << "---- model - piniocchio ------" << std::endl;
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
      std::cout << "---- model - piniocchio ------" << std::endl;
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
