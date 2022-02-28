//
// Copyright (c) 2014 CNRS
// Authors: Florent Lamiraux
//
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

#include <hpp/constraints/relative-com.hh>

#include <boost/serialization/vector.hpp>

#include <pinocchio/serialization/eigen.hpp>

#include <hpp/util/indent.hh>
#include <hpp/util/debug.hh>
#include <hpp/util/serialization.hh>

#include <hpp/pinocchio/util.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/center-of-mass-computation.hh>
#include <hpp/pinocchio/liegroup-element.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/serialization.hh>

#include <hpp/constraints/macros.hh>

namespace hpp {
  namespace constraints {

    namespace {
      static size_type size (std::vector<bool> mask)
      {
        size_type res = 0;
        for (std::vector<bool>::iterator it = mask.begin ();
            it != mask.end (); ++it)
          if (*it) ++res;
        return res;
      }
    } // namespace

    RelativeComPtr_t RelativeCom::create (const std::string& name,
                                          const DevicePtr_t& robot,
					  const JointPtr_t& joint,
					  const vector3_t reference,
                                          std::vector <bool> mask)
    {
      CenterOfMassComputationPtr_t comc =
        CenterOfMassComputation::create (robot);
      comc->add (robot->rootJoint ());
      return create (name, robot, comc, joint, reference, mask);
    }

    RelativeComPtr_t RelativeCom::create (
        const DevicePtr_t& robot,
        const CenterOfMassComputationPtr_t& comc,
        const JointPtr_t& joint, const vector3_t reference,
        std::vector <bool> mask)
    {
      return create ("RelativeCom", robot, comc, joint, reference, mask);
    }

    RelativeComPtr_t RelativeCom::create (
        const std::string& name,
        const DevicePtr_t& robot,
        const CenterOfMassComputationPtr_t& comc,
        const JointPtr_t& joint, const vector3_t reference,
        std::vector <bool> mask)
    {
      RelativeCom* ptr = new RelativeCom (robot, comc, joint, reference, mask, name);
      RelativeComPtr_t shPtr (ptr);
      return shPtr;
    }

    RelativeCom::RelativeCom (const DevicePtr_t& robot,
        const CenterOfMassComputationPtr_t& comc,
        const JointPtr_t& joint, const vector3_t reference,
        std::vector <bool> mask,
        const std::string& name) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
                              LiegroupSpace::Rn (size (mask)), name),
      robot_ (robot), comc_ (comc), joint_ (joint), reference_ (reference),
      mask_ (mask), nominalCase_ (false), jacobian_
      (3, robot->numberDof()-robot->extraConfigSpace().dimension())
    {
      if (mask[0] && mask[1] && mask[2])
        nominalCase_ = true;
      jacobian_.setZero ();
    }

    std::ostream& RelativeCom::print (std::ostream& o) const
    {
      return o << "RelativeCom: " << name () << incindent
        << iendl << "Joint: "        << (joint_ ? joint_->name() : "World")
        << iendl << "Reference: " << one_line (reference_)
        << iendl << "mask: ";
      for (size_type i=0; i<3; ++i) o << mask_ [i] << ", ";
      return o << decindent;
    }

    void RelativeCom::impl_compute (LiegroupElementRef result,
				    ConfigurationIn_t argument)
      const
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      comc_->compute (hpp::pinocchio::COM);
      const Transform3f& M = joint_->currentTransformation ();
      const vector3_t& x = comc_->com ();
      const matrix3_t& R = M.rotation ();
      const vector3_t& t = M.translation ();

      if (nominalCase_)
        result.vector () = R.transpose() * (x - t) - reference_;
      else {
        const vector3_t res ( R.transpose() * (x - t) - reference_);
        size_t index = 0;
        for (size_t i = 0; i < 3; ++i)
          if (mask_[i]) {
            result.vector () [index] = res [i];
            index++;
          }
      }
    }

    void RelativeCom::impl_jacobian (matrixOut_t jacobian,
				     ConfigurationIn_t arg) const
    {
      robot_->currentConfiguration (arg);
      robot_->computeForwardKinematics ();
      comc_->compute (hpp::pinocchio::COMPUTE_ALL);
      const ComJacobian_t& Jcom = comc_->jacobian ();
      const JointJacobian_t& Jjoint (joint_->jacobian ());
      const Transform3f& M = joint_->currentTransformation ();
      const matrix3_t& R (M.rotation ());
      const vector3_t& x (comc_->com ());
      const vector3_t& t (M.translation ());

      // Right part
      jacobian.rightCols (jacobian.cols () - Jjoint.cols ()).setZero ();
      // Left part
      // J = 0RTj ( Jcom + [ x - 0tj ]x 0Rj jJwj - 0Rj jJtj)
      jacobian_ = R.transpose() * Jcom;
      jacobian_.noalias() += (R.transpose() * R.colwise().cross(t-x)) * Jjoint.bottomRows<3>();

      if (nominalCase_) {
        jacobian.leftCols (Jjoint.cols ()).noalias() = jacobian_ - Jjoint.topRows<3>();
      } else {
        size_t index = 0;
        for (size_t i = 0; i < 3; ++i)
          if (mask_[i]) {
            jacobian.row(index).head(Jjoint.cols()) = jacobian_.row (i) - Jjoint.row(i);
            index++;
          }
      }
      hppDnum (info, "Jcom = " << std::endl << Jcom);
      hppDnum (info, "Jw = " << std::endl << Jjoint.bottomRows<3>());
      hppDnum (info, "Jv = " << std::endl << Jjoint.topRows<3>());
    }

    template<class Archive>
    void RelativeCom::serialize(Archive & ar, const unsigned int version)
    {
      using namespace boost::serialization;
      (void) version;
      ar & make_nvp("base", base_object<DifferentiableFunction>(*this));
      ar & BOOST_SERIALIZATION_NVP(robot_);
      ar & BOOST_SERIALIZATION_NVP(comc_);
      ar & BOOST_SERIALIZATION_NVP(joint_);
      ar & BOOST_SERIALIZATION_NVP(reference_);
      ar & BOOST_SERIALIZATION_NVP(mask_);
      if (!Archive::is_saving::value)
        nominalCase_ = (mask_[0] && mask_[1] && mask_[2]);
    }

    HPP_SERIALIZATION_IMPLEMENT(RelativeCom);
  } // namespace constraints
} // namespace hpp

BOOST_CLASS_EXPORT(hpp::constraints::RelativeCom)
