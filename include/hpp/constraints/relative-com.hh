//
// Copyright (c) 2014 CNRS
// Authors: Florent Lamiraux
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

#ifndef HPP_CONSTRAINTS_RELATIVE_COM_HH
# define HPP_CONSTRAINTS_RELATIVE_COM_HH

# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/config.hh>
# include <hpp/constraints/fwd.hh>

namespace hpp {
  namespace eigen {
    typedef Eigen::Matrix <double, 3, 3> matrix3_t;
    typedef Eigen::Matrix <double, 3, 1> vector3_t;
  } // namespace eigen
  namespace constraints {
    /// \addtogroup constraints
    /// \{

    /// Constraint on the relative position of the center of mass
    ///
    /// The value of the function is defined as the position of the center
    /// of mass in the reference frame of a joint.
    /**
     *  \f{eqnarray*}
     *  \mathbf{f}(\mathbf{q}) &=& R^T \left(\mathbf{x} - \mathbf{t}\right)
     *  - \mathbf{x}^{*}\\
     *  \mathbf{\dot{f}} &=& R^T \left(
     *  J_{com} + [\mathbf{x}-\mathbf{t}]_{\times}J_{joint}^{\omega}
     *  - J_{joint}^{\mathbf{v}}\right)\mathbf{\dot{q}}
     *  \f}
    **/
    /// where
    /// \li \f$
    ///     \left(\begin{array}{cc} R & \mathbf{t} \\ 0 & 1\end{array}\right)
    ///     \f$
    /// is the position of the joint,
    /// \li \f$\mathbf{x}\f$ is the position of the center of mass,
    /// \li \f$\mathbf{x}^{*}\f$ is the desired position of the center of mass
    ///     expressed in joint frame.
    class HPP_CONSTRAINTS_DLLAPI RelativeCom : public DifferentiableFunction
    {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      /// Return a shared pointer to a new instance
      static RelativeComPtr_t create (const DevicePtr_t& robot,
				      const JointPtr_t& joint,
				      const vector3_t reference);
      virtual ~RelativeCom () throw () {}
      RelativeCom (const DevicePtr_t& robot, const JointPtr_t& joint,
		   const vector3_t reference);
    protected:
      /// Compute value of error
      ///
      /// \param argument configuration of the robot,
      /// \retval result error vector
      virtual void impl_compute	(vectorOut_t result,
				 ConfigurationIn_t argument)
	const throw ();
      virtual void impl_jacobian (matrixOut_t jacobian,
				  ConfigurationIn_t arg) const throw ();
    private:
      DevicePtr_t robot_;
      JointPtr_t joint_;
      vector3_t reference_;
      mutable eigen::matrix3_t cross_;
    }; // class RelativeCom
    /// \}
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_RELATIVE_COM_HH
