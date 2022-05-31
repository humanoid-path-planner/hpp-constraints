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

#ifndef HPP_CONSTRAINTS_CONFIGURATION_CONSTRAINT_HH
#define HPP_CONSTRAINTS_CONFIGURATION_CONSTRAINT_HH

#include <Eigen/Core>
#include <hpp/constraints/config.hh>
#include <hpp/constraints/differentiable-function.hh>
#include <hpp/constraints/fwd.hh>

namespace hpp {
namespace constraints {

/// Square distance between input configuration and reference configuration
class HPP_CONSTRAINTS_DLLAPI ConfigurationConstraint
    : public DifferentiableFunction {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /// Return a shared pointer to a new instance
  static ConfigurationConstraintPtr_t create(
      const std::string& name, const DevicePtr_t& robot, ConfigurationIn_t goal,
      std::vector<bool> mask = std::vector<bool>(0));

  /// Return a shared pointer to a new instance
  /// \param weights vector of size robot->numberDof()
  static ConfigurationConstraintPtr_t create(const std::string& name,
                                             const DevicePtr_t& robot,
                                             ConfigurationIn_t goal,
                                             const vector_t& weights);

  virtual ~ConfigurationConstraint() {}

  /// \param weights vector of size robot->numberDof()
  ConfigurationConstraint(const std::string& name, const DevicePtr_t& robot,
                          ConfigurationIn_t goal, const vector_t& weights);

  const vector_t& weights() const { return weights_; }

  void weights(const vector_t& ws);

  const LiegroupElement& goal() const { return goal_; }

 protected:
  /// Compute value of error
  ///
  /// \param argument configuration of the robot,
  /// \retval result error vector
  virtual void impl_compute(LiegroupElementRef result,
                            ConfigurationIn_t argument) const;

  virtual void impl_jacobian(matrixOut_t jacobian, ConfigurationIn_t arg) const;

  std::ostream& print(std::ostream& o) const;

  bool isEqual(const DifferentiableFunction& other) const {
    const ConfigurationConstraint& castother =
        dynamic_cast<const ConfigurationConstraint&>(other);
    if (!DifferentiableFunction::isEqual(other)) return false;

    if (robot_ != castother.robot_) return false;
    if (goal_.vector() != castother.goal_.vector()) return false;
    if (weights_ != castother.weights_) return false;

    return true;
  }

 private:
  typedef Eigen::Array<bool, Eigen::Dynamic, 1> EigenBoolVector_t;
  DevicePtr_t robot_;
  LiegroupElement goal_;
  vector_t weights_;
};  // class ConfigurationConstraint
}  // namespace constraints
}  // namespace hpp
#endif  // HPP_CONSTRAINTS_CONFIGURATION_CONSTRAINT_HH
