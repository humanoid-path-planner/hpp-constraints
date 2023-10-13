// Copyright (c) 2015 - 2018, LAAS-CNRS
// Authors: Florent Lamiraux, Joseph Mirabel
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

#ifndef HPP_CONSTRAINTS_LOCKED_JOINT_HH
#define HPP_CONSTRAINTS_LOCKED_JOINT_HH

#include <hpp/constraints/explicit.hh>
#include <hpp/pinocchio/joint.hh>

namespace hpp {
namespace constraints {

/// \addtogroup constraints
/// \{

/**
Implementation of constraint specific to a locked joint.

The implicit formulation as defined in class Implicit is given by
\f{equation}
h (\mathbf{q}) = q_{out} - value
\f}
where \f$value\f$ is an element of the configuration space of the locked
joint passed to method \link
LockedJoint::create(const JointPtr_t& joint, const LiegroupElement& value)
create\endlink.

Note that \f$h\f$ takes values in \f$\mathbf{R}^{nv}\f$ where
\f$nv\f$ is the dimension of the joint tangent space.

The explicit formulation is given by
\f{equation}
q_{out} = value + rhs
\f}
where coordinates of \f$rhs\f$ corresponding to comparison types different
from Equality are set to 0.

As such, the relation between the explicit formulation and the implicit
formulation is the default one.
*/
class HPP_CONSTRAINTS_DLLAPI LockedJoint : public Explicit {
 public:
  /// Copy object and return shared pointer to copy
  virtual ImplicitPtr_t copy() const;

  /// Create locked joint and return shared pointer
  /// \param joint joint that is locked,
  /// \param value of the constant joint config,
  static LockedJointPtr_t create(const JointPtr_t& joint,
                                 const LiegroupElement& value);

  /// Create partial locked joint (only some degrees of freedom)
  /// \param joint joint that is locked,
  /// \param index first locked degree of freedom in the joint,
  /// \param value of the constant joint partial config, size of value
  ///        determines the number of degrees of freedom locked.
  /// \note valid only for translation joints.
  static LockedJointPtr_t create(const JointPtr_t& joint, const size_type index,
                                 vectorIn_t value);

  /// Create locked degrees of freedom of extra config space
  /// \param dev robot
  /// \param index index of the first  locked extra degree of freedom,
  /// \param value of the locked degrees of freedom, size of value
  ///        determines the number of degrees of freedom locked.
  static LockedJointPtr_t create(const DevicePtr_t& dev, const size_type index,
                                 vectorIn_t value);

  /// Return shared pointer to copy
  /// \param other instance to copy.
  static LockedJointPtr_t createCopy(LockedJointConstPtr_t other);

  /// Get index of locked degree of freedom in robot configuration vector
  size_type rankInConfiguration() const;

  /// Get index of locked degree of freedom in robot velocity vector.
  size_type rankInVelocity() const;

  /// Get the configuration size of the joint.
  size_type configSize() const;

  /// Get number of degrees of freedom of the joint
  size_type numberDof() const;

  /// Get configuration space of locked joint
  const LiegroupSpacePtr_t& configSpace() const;

  /// Get the value of the locked joint.
  vectorIn_t value() const;

  /// Set the value of the locked joint.
  void value(vectorIn_t value);

  /// Return shared pointer to joint
  const JointPtr_t& joint() { return joint_; }
  /// Return the joint name.
  const std::string& jointName() const { return jointName_; }
  /// Print object in a stream
  std::ostream& print(std::ostream& os) const;

  /// Get pair of joints whose relative pose is fully constrained
  ///
  /// \param robot the device that this constraint is applied to
  /// \return pair of pointers to the parent joint and the locked joint,
  /// arranged in order of increasing joint index
  /// \note absolute pose is considered relative pose with respect to
  /// "universe". "universe" is returned as a nullpointer
  /// as the first element of the pair, if applicable.
  virtual std::pair<JointConstPtr_t, JointConstPtr_t>
  doesConstrainRelPoseBetween(DeviceConstPtr_t robot) const;

 protected:
  /// Constructor
  /// \param joint joint that is locked,
  /// \param value of the constant joint config,
  LockedJoint(const JointPtr_t& joint, const LiegroupElement& value);
  /// Constructor of partial locked joint
  /// \param joint joint that is locked,
  /// \param index first locked degree of freedom in the joint,
  /// \param value of the constant joint partial config, size of value
  ///        determines the number of degrees of freedom locked.
  /// \note valid only for translation joints.
  LockedJoint(const JointPtr_t& joint, const size_type index, vectorIn_t value);
  /// Constructor of locked degrees of freedom of extra config space
  /// \param robot robot
  /// \param index index of the first  locked extra degree of freedom,
  /// \param value of the locked degrees of freedom, size of value
  ///        determines the number of degrees of freedom locked.
  LockedJoint(const DevicePtr_t& robot, const size_type index,
              vectorIn_t value);
  /// Copy constructor
  LockedJoint(const LockedJoint& other);
  /// Test equality with other instance
  /// \param other object to copy
  /// \param swapAndTest whether we should also check other == this
  virtual bool isEqual(const Implicit& other, bool swapAndTest) const;

  void init(const LockedJointPtr_t& self);

 private:
  std::string jointName_;
  JointPtr_t joint_;
  LiegroupSpacePtr_t configSpace_;
  /// Weak pointer to itself
  LockedJointWkPtr_t weak_;

  LockedJoint() {}
  HPP_SERIALIZABLE();
};  // class LockedJoint

/// \}
inline std::ostream& operator<<(std::ostream& os, const LockedJoint& lj) {
  return lj.print(os);
}
}  // namespace constraints
}  // namespace hpp

BOOST_CLASS_EXPORT_KEY(hpp::constraints::LockedJoint)
#endif  // HPP_CONSTRAINTS_LOCKED_JOINT_HH
