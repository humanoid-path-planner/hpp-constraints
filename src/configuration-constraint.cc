// Copyright (c) 2015, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
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

#include <hpp/constraints/configuration-constraint.hh>

#include <pinocchio/multibody/liegroup/liegroup.hpp>
#include <pinocchio/multibody/model.hpp>

#include <hpp/util/debug.hh>
#include <hpp/util/indent.hh>
#include <hpp/pinocchio/util.hh>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/liegroup-element.hh>
#include <hpp/pinocchio/joint-collection.hh>

namespace hpp {
  namespace constraints {

    ConfigurationConstraintPtr_t ConfigurationConstraint::create (
        const std::string& name, const DevicePtr_t& robot,
        ConfigurationIn_t goal, std::vector <bool> mask)
    {
      vector_t ws (vector_t::Ones(robot->numberDof ()));
      for (std::size_t i = 0; i < mask.size (); ++i) {
        if (!mask[i]) ws[i] = 0;
      }

      ConfigurationConstraint* ptr = new ConfigurationConstraint
        (name, robot, goal, ws);
      return ConfigurationConstraintPtr_t (ptr);
    }

    ConfigurationConstraintPtr_t ConfigurationConstraint::create (
        const std::string& name, const DevicePtr_t& robot,
        ConfigurationIn_t goal, const vector_t& weights)
    {
      ConfigurationConstraint* ptr = new ConfigurationConstraint
        (name, robot, goal, weights);
      return ConfigurationConstraintPtr_t (ptr);
    }

    ConfigurationConstraint::ConfigurationConstraint (
        const std::string& name, const DevicePtr_t& robot,
        ConfigurationIn_t goal, const vector_t& ws) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
                              LiegroupSpace::R1 (), name),
      robot_ (robot), weights_ ()
    {
      weights(ws);
      LiegroupSpacePtr_t s (LiegroupSpace::createCopy(robot->configSpace()));
      s->mergeVectorSpaces();
      goal_ = LiegroupElement (goal, s);
    }

    void ConfigurationConstraint::weights (const vector_t& ws)
    {
      if (ws.size() != robot_->numberDof())
        throw std::invalid_argument("Size of weights vector should be the same "
            "as the robot number DoFs.");
      weights_ = ws;

      activeParameters_ = ArrayXb::Constant (inputSize(), true);
      activeDerivativeParameters_ = ArrayXb::Constant (weights_.size(), true);

      const pinocchio::Model& model = robot_->model();
      for (int i = 1; i < model.njoints; ++i) {

        if ((weights_.segment(model.joints[i].idx_v(), model.joints[i].nv())
              .array() == 0).all()) {
          activeParameters_.segment(
              model.joints[i].idx_q(), model.joints[i].nq()).setConstant(false);
          activeDerivativeParameters_.segment(
              model.joints[i].idx_v(), model.joints[i].nv()).setConstant(false);
        }
      }
    }

    std::ostream& ConfigurationConstraint::print (std::ostream& o) const
    {
      o << "ConfigurationConstraint: " << name() << incindent << iendl
        << "weights: " << one_line (weights_);
      return o << decindent;
    }

    void ConfigurationConstraint::impl_compute (LiegroupElementRef result,
                                                ConfigurationIn_t argument)
      const
    {
      using namespace hpp::pinocchio;
      LiegroupElementConstRef a (argument, goal_.space());
      result.vector () [0] = 0.5 * weights_.dot((goal_ - a).cwiseAbs2());
    }

    void ConfigurationConstraint::impl_jacobian (matrixOut_t jacobian,
        ConfigurationIn_t argument) const
    {
      using namespace hpp::pinocchio;

      LiegroupElementConstRef a (argument, goal_.space());
      jacobian.leftCols (robot_->numberDof ()).noalias() = (goal_ - a).transpose();

      // Apply jacobian of the difference on the right.
      goal_.space()->dDifference_dq0<pinocchio::InputTimesDerivative>
        (argument, goal_.vector(),
          jacobian.leftCols (robot_->numberDof ()));

      jacobian.leftCols (robot_->numberDof ()).array()
        *= weights_.array().transpose();
    }
  } // namespace constraints
} // namespace hpp
