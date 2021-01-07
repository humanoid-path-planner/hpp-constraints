// Copyright (c) 2015, LAAS-CNRS
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//
// This file is part of hpp-core.
// hpp-core is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-core is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-core. If not, see <http://www.gnu.org/licenses/>.

#include "hpp/constraints/locked-joint.hh"

#include <sstream>

#include <boost/assign/list_of.hpp>
#include <boost/serialization/weak_ptr.hpp>

#include <hpp/util/debug.hh>
#include <hpp/util/serialization.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>

#include <hpp/constraints/affine-function.hh>
#include <hpp/constraints/explicit/implicit-function.hh>

namespace hpp {
  namespace constraints {
    using boost::assign::list_of;
    namespace {
        template <typename T>
	std::string numToStr (const T& v) {
 	  std::stringstream ss; ss<<v; return ss.str ();
	}
    }

    /// Copy object and return shared pointer to copy
    ImplicitPtr_t LockedJoint::copy () const
    {
      return createCopy (weak_.lock ());
    }

    LockedJointPtr_t LockedJoint::create (const JointPtr_t& joint,
					  const LiegroupElement& value)
    {
      LockedJoint* ptr = new LockedJoint (joint, value);
      LockedJointPtr_t shPtr (ptr);
      ptr->init (shPtr);
      return shPtr;
    }

    LockedJointPtr_t LockedJoint::create (const JointPtr_t& joint,
                                          const size_type index,
					  vectorIn_t value)
    {
      LockedJoint* ptr = new LockedJoint (joint, index, value);
      LockedJointPtr_t shPtr (ptr);
      ptr->init (shPtr);
      return shPtr;
    }

    LockedJointPtr_t LockedJoint::create (const DevicePtr_t& dev,
        const size_type index, vectorIn_t value)
    {
      LockedJoint* ptr = new LockedJoint (dev, index, value);
      LockedJointPtr_t shPtr (ptr);
      ptr->init (shPtr);
      return shPtr;
    }

    LockedJointPtr_t LockedJoint::createCopy (LockedJointConstPtr_t other)
    {
      LockedJoint* ptr = new LockedJoint (*other);
      LockedJointPtr_t shPtr (ptr);
      ptr->init (shPtr);
      return shPtr;
    }

    size_type LockedJoint::rankInConfiguration () const
    {
      return outputConf_[0].first;
    }

    size_type LockedJoint::rankInVelocity () const
    {
      return outputVelocity_[0].first;
    }

    size_type LockedJoint::configSize () const
    {
      return configSpace_->nq ();
    }

    size_type LockedJoint::numberDof () const
    {
      return configSpace_->nv ();
    }

    const LiegroupSpacePtr_t& LockedJoint::configSpace () const
    {
      return configSpace_;
    }

    void LockedJoint::rightHandSideFromConfig (ConfigurationIn_t config)
    {
      Implicit::rightHandSideFromConfig (config);
#ifndef NDEBUG
      vector_t rhs (rightHandSide ());
      if (parameterSize () != 0) {
        LiegroupElement q (config.segment (rankInConfiguration(), configSize ()),
                           configSpace_);
        rightHandSide (q - configSpace_->neutral ());
      }
      assert (rhs == rightHandSide ());
#endif
    }

    LockedJoint::LockedJoint (const JointPtr_t& joint,
                              const LiegroupElement& value) :
      Explicit (joint->robot()->configSpace (),
                ConstantFunction::create(value, 0, 0, "LockedJoint " +
                                         joint->name()+' '+
                                         numToStr(value.vector().transpose())),
                segments_t(), // input conf
                list_of(segment_t (joint->rankInConfiguration(),
                                   joint->configSize())), // output conf
                segments_t(), // input vel
                list_of(segment_t (joint->rankInVelocity     (),
                                   joint->numberDof ())), // output vel
                ComparisonTypes_t(joint->numberDof(), Equality),
		std::vector<bool>(joint->numberDof(), true)),
      jointName_ (joint->name ()),
      configSpace_ (joint->configurationSpace ())
    {
      assert (HPP_DYNAMIC_PTR_CAST (explicit_::ImplicitFunction,
                                    functionPtr ()));
      assert (rightHandSideSize () == joint->numberDof ());
      assert (*(value.space ()) == *configSpace_);
    }

    LockedJoint::LockedJoint (const JointPtr_t& joint, const size_type index,
                              vectorIn_t value) :
      Explicit (joint->robot()->configSpace (),
                ConstantFunction::create
                (LiegroupElement(value, LiegroupSpace::Rn(joint->configSize () -
                                                          index)), 0, 0,
                 "partial " + joint->name()),
                segments_t(), // input conf
                segments_t(), // input vel
                list_of(segment_t (joint->rankInConfiguration(),
                                   joint->configSize()-index)), // output conf
                list_of(segment_t (joint->rankInVelocity     (),
                                   joint->numberDof ()-index)), // output vel
                ComparisonTypes_t(joint->numberDof()-index, Equality),
		std::vector<bool>(joint->numberDof(), true)),
      jointName_ ("partial_" + joint->name ()),
      joint_ (joint),
      configSpace_ (LiegroupSpace::Rn (joint->configSize () - index))
    {
      assert (HPP_DYNAMIC_PTR_CAST (explicit_::ImplicitFunction,
                                    functionPtr ()));
      assert (joint->numberDof () == joint->configSize ());
      // rightHandSide (value);
      assert (rightHandSideSize () == value.size());
    }

    LockedJoint::LockedJoint (const DevicePtr_t& dev, const size_type index,
                              vectorIn_t value) :
      Explicit (dev->configSpace (),
                ConstantFunction::create
                (LiegroupElement(value, LiegroupSpace::Rn (value.size())), 0, 0,
                 dev->name() + "_extraDof" + numToStr (index)),
                segments_t(), // input conf
                segments_t(), // input vel
                list_of(segment_t (dev->configSize () -
                                   dev->extraConfigSpace().dimension() + index,
                                   value.size())), // output conf
                list_of(segment_t (dev->numberDof ()  -
                                   dev->extraConfigSpace().dimension() + index,
                                   value.size())), // output vel
                ComparisonTypes_t(value.size(), Equality),
		std::vector<bool>(value.size(), true)),
      jointName_ (dev->name() + "_extraDof" + numToStr (index)),
      joint_ (JointPtr_t ()), configSpace_ (LiegroupSpace::Rn (value.size ()))
    {
      assert (HPP_DYNAMIC_PTR_CAST (explicit_::ImplicitFunction,
                                    functionPtr ()));
      assert (value.size() > 0);
      assert (rankInConfiguration() + value.size() <= dev->configSize());
      // rightHandSide (value);
      assert (rightHandSideSize () == value.size());
    }

    void LockedJoint::init (const LockedJointPtr_t& self)
    {
      Explicit::init(self);
      weak_ = self;
    }

    std::ostream& LockedJoint::print (std::ostream& os) const
    {
      LiegroupElement v; vector_t empty;
      function().value(v, empty);
      os << "Locked joint " << jointName_
	 << ", value = " << pinocchio::displayConfig (v.vector())
        << ": rank in configuration = " << rankInConfiguration()
        << ": rank in velocity = " << rankInVelocity()
        << std::endl;
      return os;
    }

    LockedJoint::LockedJoint (const LockedJoint& other) :
      Explicit (other), jointName_ (other.jointName_),
      joint_ (other.joint_), configSpace_ (other.configSpace_), weak_ ()
    {
    }

    bool LockedJoint::isEqual (const Implicit& other, bool swapAndTest) const
    {
      try {
	const LockedJoint& lj =
	  dynamic_cast <const LockedJoint&> (other);
	if (!Implicit::isEqual (other, false)) return false;
	if (jointName_ != lj.jointName_) return false;
	if (rankInConfiguration() != lj.rankInConfiguration()) return false;
	if (rankInVelocity() != lj.rankInVelocity()) return false;
	if (*configSpace_ != *(lj.configSpace_)) return false;
	if (swapAndTest) return lj.isEqual (*this, false);
	return true;
      } catch (const std::bad_cast& err) {
	return false;
      }
    }

    template<class Archive>
    void LockedJoint::serialize(Archive & ar, const unsigned int version)
    {
      (void) version;
      ar & boost::serialization::make_nvp("base",
           boost::serialization::base_object<Explicit>(*this));
      ar & BOOST_SERIALIZATION_NVP(jointName_);
      ar & BOOST_SERIALIZATION_NVP(joint_);
      ar & BOOST_SERIALIZATION_NVP(configSpace_);
      ar & BOOST_SERIALIZATION_NVP(weak_);
    }

    HPP_SERIALIZATION_IMPLEMENT(LockedJoint);
  } // namespace constraints
} // namespace hpp

BOOST_CLASS_EXPORT_IMPLEMENT(hpp::constraints::LockedJoint)
