// Copyright (c) 2018 CNRS
// Authors: Joseph Mirabel
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

#ifndef HPP_CONSTRAINTS_MANIPULABILITY_HH
#define HPP_CONSTRAINTS_MANIPULABILITY_HH

#include <hpp/constraints/config.hh>
#include <hpp/constraints/differentiable-function.hh>
#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/matrix-view.hh>

namespace hpp {
namespace constraints {
HPP_PREDEF_CLASS(Manipulability);
typedef shared_ptr<Manipulability> ManipulabilityPtr_t;

/// \addtogroup constraints
/// \{

/// Manipulability index
///
/// This function takes as input another differentiable function \f$f_1\f$
/// defined over the same configuration space and computes a one dimensional
/// value as follows:
///
/// \f[
/// f(\mathbf{q}) = \log\det(J_1 J_1^T)
/// \f]
/// where \f$J_1\f$ is the Jacobian matrix of \f$f_1\f$.
///
/// Inserting a \link hpp::constraints::Implicit constraint  \endlink with this
/// function with comparison type EQUALITY into a numerical solver will
/// make the solution reach a manipulability index specified by the right hand
/// side of the constraint.
///
/// \par Lock joints
///
/// In some cases, it is useful not to consider some joints in the kinematic
/// chain. For example, when evaluating the manipulability of a robotic arm
/// moving on a prismatic rail, it can be useful to consider the manipulability
/// of the system when the rail is locked. To do so, call method
/// Manipulability::lockJoint.
///
/// \note The Jacobian of this function is computed by finite difference.

class HPP_CONSTRAINTS_DLLAPI Manipulability : public DifferentiableFunction {
 public:
  virtual ~Manipulability() {}

  static ManipulabilityPtr_t create(DifferentiableFunctionPtr_t function,
                                    DevicePtr_t robot, std::string name) {
    return ManipulabilityPtr_t(new Manipulability(function, robot, name));
  }

  /// Lock a joint
  ///
  /// Consider this joint as fixed when computing the manipulability
  void lockJoint(const JointPtr_t& joint);

  /// Get robot
  const DevicePtr_t& robot() const { return robot_; }

 protected:
  /// \brief Concrete class constructor should call this constructor.
  ///
  /// \param function the function which must be analysed
  /// \param name function's name
  Manipulability(DifferentiableFunctionPtr_t function, DevicePtr_t robot,
                 std::string name);

  void impl_compute(LiegroupElementRef result, vectorIn_t argument) const;

  void impl_jacobian(matrixOut_t jacobian, vectorIn_t arg) const;

  bool isEqual(const DifferentiableFunction& other) const {
    const Manipulability& castother =
        dynamic_cast<const Manipulability&>(other);
    if (!DifferentiableFunction::isEqual(other)) return false;

    if (function_ != castother.function_) return false;
    if (robot_ != castother.robot_) return false;
    if (cols_.cols() != castother.cols_.cols()) return false;
    if (J_ != castother.J_) return false;
    if (J_JT_ != castother.J_JT_) return false;

    return true;
  }

 private:
  DifferentiableFunctionPtr_t function_;
  DevicePtr_t robot_;

  Eigen::ColBlockIndices cols_;
  ArrayXb activeAndNonLockedDerivParams_;
  mutable matrix_t J_, J_JT_;
};  // class Manipulability
/// \}
}  // namespace constraints
}  // namespace hpp

#endif  // HPP_CONSTRAINTS_MANIPULABILITY_HH
