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

#ifndef HPP_CONSTRAINTS_RELATIVE_POSITION_HH
# define HPP_CONSTRAINTS_RELATIVE_POSITION_HH

# include <boost/assign/list_of.hpp>
# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/config.hh>

namespace hpp {
  namespace eigen {
    typedef Eigen::Matrix <double, 3, 3> matrix3_t;
    typedef Eigen::Matrix <double, 3, 1> vector3_t;
  } // namespace eigen
  namespace constraints {
    /// \addtogroup constraints
    /// \{

    /// Constraint on the relative position of two robot joints
    ///
    /// The value of the function is defined as follows:
    /** \f{eqnarray*}
     *
     *  \mathbf{r} (\mathbf{q}) &=& \mathbf{x}_1 - M_1^{-1}(\mathbf{q})
     *  M_2 (\mathbf{q})\mathbf{x}_2\\
     *  &=& \mathbf{x}_1
     *  + R_1^T (\mathbf{q}) (\mathbf{t}_1 (\mathbf{q}) - \mathbf{t}_2 (\mathbf{q}))
     *  - R_1^T (\mathbf{q}) R_2 (\mathbf{q}) \mathbf{x}_2 \\
     *  \mathbf{\dot{e}} &=& R_1^T (\mathbf{q}) \left(
     *    \left[R_2 (\mathbf{q})\mathbf{x}_2\right]_{\times} (\omega_2 - \omega_1)
     *  + \left[\mathbf{t}_1 - \mathbf{t}_2\right]_{\times} \omega_1
     *  + \mathbf{\dot{t}}_1 - \mathbf{\dot{t}}_2 \right)\\
     *  &=& R_1^T (\mathbf{q}) \left(
     *    \left[R_2 (\mathbf{q})\mathbf{x}_2\right]_{\times} \omega_2
     *  + \left[\mathbf{t}_1 - (R_2 (\mathbf{q})\mathbf{x}_2 + \mathbf{t}_2)\right]_{\times} \omega_1
     *  + \mathbf{\dot{t}}_1 - \mathbf{\dot{t}}_2 \right)\\
     *  &=& R_1^T (\mathbf{q}) \left(
     *    \left[R_2 (\mathbf{q})\mathbf{x}_2\right]_{\times} J_{joint\ 2}^{\omega}(\mathbf{q})
     *  + \left[\mathbf{t}_1 - (R_2 (\mathbf{q})\mathbf{x}_2 + \mathbf{t}_2)\right]_{\times} J_{joint\ 1}^{\omega}(\mathbf{q})
     *  + J_{joint\ 1}^{\mathbf{v}}(\mathbf{q}) - J_{joint\ 2}^{\mathbf{v}}(\mathbf{q}) \right)
     *  \mathbf{\dot{q}}\\
     *  \f}
    **/
    /// where \f$\mathbf{x}_1\f$ and \f$\mathbf{x}_2\f$ are the position of
    /// the points to match expressed in the local frame of each joint.
    class HPP_CONSTRAINTS_DLLAPI RelativePosition :
      public DifferentiableFunction
    {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      /// Return a shared pointer to a new instance
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      /// \param joint1, joint2 joints the relative positions of which is
      ///                       constrained.
      /// \param point1, point2 the points on each joint that should coincide
      ///   expressed in joint local frames.
      /// \param mask which component of the error vector to take into
      ///        account.
      static RelativePositionPtr_t create
	(const std::string& name, const DevicePtr_t& robot,
         const JointPtr_t& joint1, const JointPtr_t& joint2,
         const vector3_t& point1, const vector3_t& point2,
	 std::vector <bool> mask = boost::assign::list_of (true)(true)(true));

      /// Return a shared pointer to a new instance
      ///
      /// \param robot the robot the constraints is applied to,
      /// \param joint1, joint2 joints the relative positions of which is
      ///                       constrained.
      /// \param point1, point2 the points on each joint that should coincide
      ///   expressed in joint local frames.
      /// \param mask which component of the error vector to take into
      ///        account.
      static RelativePositionPtr_t create
	(const DevicePtr_t& robot, const JointPtr_t& joint1,
	 const JointPtr_t& joint2, const vector3_t& point1,
	 const vector3_t& point2,
	 std::vector <bool> mask = boost::assign::list_of (true)(true)(true));
      virtual ~RelativePosition () throw () {}
      /// Constructor
      ///
      /// \param robot the robot the constraints is applied to,
      /// \param joint1, joint2 joints the relative positions of which is
      ///                       constrained.
      /// \param point1, point2 the points on each joint that should coincide
      ///   expressed in joint local frames.
      /// \param mask which component of the error vector to take into
      ///        account.
      RelativePosition (const std::string& name, const DevicePtr_t& robot,
                        const JointPtr_t& joint1, const JointPtr_t& joint2,
                        const vector3_t& point1, const vector3_t& point2,
                        std::vector <bool> mask = boost::assign::list_of (true)(true)(true));

      /// Get reference point in joint 1
      const vector3_t& pointInJoint1 () const
      {
	return point1_;
      }
      /// Set reference point in joint 1
      void pointInJoint1 (const vector3_t& reference)
      {
	point1_ = reference;
      }
      /// Get reference point in joint 2
      const vector3_t& pointInJoint2 () const
      {
	return point2_;
      }
      /// Set reference point in joint 2
      void pointInJoint2 (const vector3_t& reference)
      {
	point2_ = reference;
      }
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
      JointPtr_t joint1_;
      JointPtr_t joint2_;
      vector3_t point1_;
      vector3_t point2_;
      std::vector <bool> mask_;
      mutable vector3_t global2_;
      mutable vector3_t point2in1_;
      mutable vector3_t global2minusT1_;
      mutable vector3_t R2x2_;
      mutable vector_t result_;
      mutable ComJacobian_t jacobian_;
      mutable eigen::matrix3_t cross1_;
      mutable eigen::matrix3_t cross2_;
      mutable eigen::matrix3_t R1T_;
      mutable eigen::vector3_t T1_;
    }; // class RelativePosition
    /// \}
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_RELATIVE_POSITION_HH
