//
// Copyright (c) 2015 CNRS
// Authors: Joseph Mirabel
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

#ifndef HPP_CONSTRAINTS_COM_BETWEEN_FEET_HH
#define HPP_CONSTRAINTS_COM_BETWEEN_FEET_HH

#include <hpp/constraints/config.hh>
#include <hpp/constraints/differentiable-function.hh>
#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/symbolic-calculus.hh>
#include <hpp/constraints/tools.hh>

namespace hpp {
namespace constraints {

/**
 *  Constraint on the relative position of the center of mass
 *
 *  The value of the function is defined as the position of the center
 *  of mass in the reference frame of a joint.
 *
 *  \f{eqnarray*}
 *  \mathbf{f}(\mathbf{q}) &=&
 *  \left(\begin{array}{c}
 *    ( x_{com} - x_{ref} ) \cdot u_z \\
 *    ( R^T (e \wedge u) ) \cdot u_z \\
 *    ( x_{com} - x_L ) \cdot (u)\\
 *    ( x_{com} - x_R ) \cdot (u)\\
 *  \end{array}\right)
 *  \f}
 *  where
 *  \li \f$\mathbf{x}_{com}\f$ is the position of the center of mass,
 *  \li \f$\mathbf{x_L}\f$ is the position of the left joint,
 *  \li \f$\mathbf{x_R}\f$ is the position of the right joint,
 *  \li \f$\mathbf{x}_{ref}\f$ is the desired position of the center of mass
 *      expressed in reference joint frame.
 *  \li \f$ u = x_R - x_L \f$
 *  \li \f$ e = x_{com} - (\frac{x_L + x_R}{2})\f$
 **/
class HPP_CONSTRAINTS_DLLAPI ComBetweenFeet : public DifferentiableFunction {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /// Return a shared pointer to a new instance
  static ComBetweenFeetPtr_t create(
      const std::string& name, const DevicePtr_t& robot,
      const JointPtr_t& jointLeft, const JointPtr_t& jointRight,
      const vector3_t pointLeft, const vector3_t pointRight,
      const JointPtr_t& jointReference, const vector3_t pointRef,
      std::vector<bool> mask = {true, true, true, true});

  /// Return a shared pointer to a new instance
  static ComBetweenFeetPtr_t create(
      const std::string& name, const DevicePtr_t& robot,
      const CenterOfMassComputationPtr_t& comc, const JointPtr_t& jointLeft,
      const JointPtr_t& jointRight, const vector3_t pointLeft,
      const vector3_t pointRight, const JointPtr_t& jointReference,
      const vector3_t pointRef,
      std::vector<bool> mask = {true, true, true, true});

  virtual ~ComBetweenFeet() {}

  ComBetweenFeet(const std::string& name, const DevicePtr_t& robot,
                 const CenterOfMassComputationPtr_t& comc,
                 const JointPtr_t& jointLeft, const JointPtr_t& jointRight,
                 const vector3_t pointLeft, const vector3_t pointRight,
                 const JointPtr_t& jointReference, const vector3_t pointRef,
                 std::vector<bool> mask);

 protected:
  /// Compute value of error
  ///
  /// \param argument configuration of the robot,
  /// \retval result error vector
  virtual void impl_compute(LiegroupElementRef result,
                            ConfigurationIn_t argument) const;

  virtual void impl_jacobian(matrixOut_t jacobian, ConfigurationIn_t arg) const;

  bool isEqual(const DifferentiableFunction& other) const {
    const ComBetweenFeet& castother =
        dynamic_cast<const ComBetweenFeet&>(other);
    if (!DifferentiableFunction::isEqual(other)) return false;

    if (robot_ != castother.robot_) return false;
    if (com_->centerOfMassComputation() !=
        castother.com_->centerOfMassComputation())
      return false;
    if (left_->joint() != castother.left_->joint()) return false;
    if (right_->joint() != castother.right_->joint()) return false;
    if (left_->local() != castother.left_->local()) return false;
    if (right_->local() != castother.right_->local()) return false;
    if (jointRef_ != castother.jointRef_) return false;
    if (pointRef_ != castother.pointRef_) return false;
    if (mask_ != castother.mask_) return false;

    return true;
  }

 private:
  DevicePtr_t robot_;
  mutable Traits<PointCom>::Ptr_t com_;
  Traits<PointInJoint>::Ptr_t left_, right_;
  eigen::vector3_t pointRef_;
  JointPtr_t jointRef_;
  typedef Difference<PointCom, PointInJoint> DiffPCPiJ;
  typedef Difference<PointInJoint, PointInJoint> DiffPiJPiJ;
  typedef Sum<PointInJoint, PointInJoint> SumPiJPiJ;
  typedef CrossProduct<Difference<PointCom, ScalarMultiply<SumPiJPiJ> >,
                       DiffPiJPiJ>
      ECrossU_t;
  mutable Traits<DiffPCPiJ>::Ptr_t xmxl_, xmxr_;
  mutable Traits<DiffPiJPiJ>::Ptr_t u_;
  mutable Traits<ECrossU_t>::Ptr_t ecrossu_;
  mutable Traits<RotationMultiply<ECrossU_t> >::Ptr_t expr_;
  mutable Traits<ScalarProduct<DiffPCPiJ, DiffPiJPiJ> >::Ptr_t xmxlDotu_,
      xmxrDotu_;
  std::vector<bool> mask_;
  mutable eigen::matrix3_t cross_;
};  // class ComBetweenFeet
}  // namespace constraints
}  // namespace hpp
#endif  // HPP_CONSTRAINTS_COM_BETWEEN_FEET_HH
