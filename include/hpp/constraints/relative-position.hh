//
// Copyright (c) 2014 CNRS
// Authors: Florent Lamiraux
//
//
// This file is part of hpp-wholebody-step.
// hpp-wholebody-step-planner is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-wholebody-step-planner is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-wholebody-step-planner. If not, see
// <http://www.gnu.org/licenses/>.

#ifndef HPP_CONSTRAINTS_RELATIVE_POSITION_HH
# define HPP_CONSTRAINTS_RELATIVE_POSITION_HH

# include <roboptim/core/differentiable-function.hh>
# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/config.hh>

namespace hpp {
  namespace constraints {
    /// Constraint on the relative position of two robot joints
    ///
    /// The value of the function is defined as follows:
    /// \f{eqnarray*}
    /// \mathbf{r} (\mathbf{q}) &=& M_1(\mathbf{q})\mathbf{x}_1 -
    /// M_2 (\mathbf{q})\mathbf{x}_2\\
    /// &=& R_1(\mathbf{q})\mathbf{x}_1 + \mathbf{t}_1 (\mathbf{q})
    /// -R_2 (\mathbf{q})\mathbf{x}_2 -  \mathbf{t}_2 (\mathbf{q})\\
    /// \mathbf{\dot{e}} &=&
    /// \left[\omega_1\right]_{\times} R_1(\mathbf{q})\mathbf{x}_1
    /// + \mathbf{\dot{t}}_1
    /// -\left[\omega_2\right]_{\times}R_2 (\mathbf{q})\mathbf{x}_2
    ///  - \mathbf{\dot{t}}_2\\
    /// &=&
    /// -\left[R_1(\mathbf{q})\mathbf{x}_1\right]_{\times} \omega_1
    /// + \mathbf{\dot{t}}_1
    /// +\left[R_2 (\mathbf{q})\mathbf{x}_2\right]_{\times}\omega_2
    ///  - \mathbf{\dot{t}}_2\\
    /// &=&
    /// \left(-\left[R_1(\mathbf{q})\mathbf{x}_1\right]_{\times}
    /// J_{joint\ 1}^{\omega}(\mathbf{q})
    /// + J_{joint\ 1}^{\mathbf{v}}(\mathbf{q})
    /// +\left[R_2 (\mathbf{q})\mathbf{x}_2\right]_{\times}J_{joint\ 2}^{\omega}
    /// (\mathbf{q}) - J_{joint\ 2}^{\mathbf{v}}(\mathbf{q})\right)
    /// \mathbf{\dot{q}}\\
    /// \f}
    /// where \f$\mathbf{x}_1\f$ and \f$\mathbf{x}_1\f$ are the position of
    /// the points to match expressed in the local frame of each joint.
    class HPP_CONSTRAINTS_DLLAPI RelativePosition :
      public DifferentiableFunction_t
    {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      ROBOPTIM_FUNCTION_FWD_TYPEDEFS (DifferentiableFunction_t);
      /// Return a shared pointer to a new instance
      ///
      /// \param robot the robot the constraints is applied to,
      /// \param joint1, joint2 joints the relative positions of which is
      ///                       constrained.
      /// \param point1, point2 the points on each joint that should coincide
      ///   expressed in joint local frames.
      static RelativePositionPtr_t create (const DevicePtr_t& robot,
					   const JointPtr_t& joint1,
					   const JointPtr_t& joint2,
					   const vector3d& point1,
					   const vector3d& point2);
      virtual ~RelativePosition () throw () {}
    protected:
      /// Constructor
      ///
      /// \param robot the robot the constraints is applied to,
      /// \param joint1, joint2 joints the relative positions of which is
      ///                       constrained.
      /// \param point1, point2 the points on each joint that should coincide
      ///   expressed in joint local frames.
      RelativePosition (const DevicePtr_t& robot,
			const JointPtr_t& joint1,
			const JointPtr_t& joint2,
			const vector3d& point1,
			const vector3d& point2);
      /// Compute value of error
      ///
      /// \param argument configuration of the robot,
      /// \retval result error vector
      virtual void impl_compute	(result_t& result,
				 const argument_t& argument) const throw ();
      virtual void impl_jacobian (jacobian_t &jacobian,
				  const argument_t &arg) const throw ();
      virtual void impl_gradient (gradient_t &gradient,
				  const argument_t &argument,
				  size_type functionId=0) const throw ();
    private:
      DevicePtr_t robot_;
      JointPtr_t joint1_;
      JointPtr_t joint2_;
      vector3d point1_;
      vector3d point2_;
      mutable vector3d global1_;
      mutable vector3d global2_;
      mutable vector3d R1x1_;
      mutable vector3d R2x2_;
      mutable matrix3d cross1_;
      mutable matrix3d cross2_;
      mutable matrix_t Jjoint1_;
      mutable matrix_t Jjoint2_;
    }; // class RelativePosition
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_RELATIVE_POSITION_HH
