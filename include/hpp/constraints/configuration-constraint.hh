//
// Copyright (c) 2015 CNRS
// Authors: Joseph Mirabel
//
//
// This file is part of hpp-constraints.
// hpp-constraints is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-constraints is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-constraints. If not, see
// <http://www.gnu.org/licenses/>.

#ifndef HPP_CONSTRAINTS_CONFIGURATION_CONSTRAINT_HH
# define HPP_CONSTRAINTS_CONFIGURATION_CONSTRAINT_HH

# include <boost/assign/list_of.hpp>
# include <Eigen/Core>
# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/config.hh>
# include <hpp/constraints/fwd.hh>

namespace hpp {
  namespace constraints {

    /// Square distance between input configuration and reference configuration
    class HPP_CONSTRAINTS_DLLAPI ConfigurationConstraint : public DifferentiableFunction
    {
      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        /// Return a shared pointer to a new instance
        static ConfigurationConstraintPtr_t create (
            const std::string& name, const DevicePtr_t& robot,
            ConfigurationIn_t goal,
            std::vector <bool> mask = std::vector <bool> (0));

        /// Return a shared pointer to a new instance
        /// \param weights vector of size robot->numberDof()
        static ConfigurationConstraintPtr_t create (
            const std::string& name, const DevicePtr_t& robot,
            ConfigurationIn_t goal,
            const vector_t& weights);

        virtual ~ConfigurationConstraint () throw () {}

        /// \param weights vector of size robot->numberDof()
        ConfigurationConstraint (const std::string& name,
            const DevicePtr_t& robot, ConfigurationIn_t goal,
            const vector_t& weights);

      protected:
        /// Compute value of error
        ///
        /// \param argument configuration of the robot,
        /// \retval result error vector
        virtual void impl_compute (LiegroupElementRef result,
                                   ConfigurationIn_t argument) const throw ();

        virtual void impl_jacobian (matrixOut_t jacobian,
            ConfigurationIn_t arg) const throw ();
      private:
        typedef Eigen::Array <bool, Eigen::Dynamic, 1> EigenBoolVector_t;
        DevicePtr_t robot_;
        LiegroupElement goal_;
        vector_t weights_;
    }; // class ComBetweenFeet
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_CONFIGURATION_CONSTRAINT_HH
