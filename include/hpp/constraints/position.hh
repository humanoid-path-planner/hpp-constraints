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

#ifndef HPP_CONSTRAINTS_POSITION_HH
# define HPP_CONSTRAINTS_POSITION_HH

# include <boost/assign/list_of.hpp>
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

    /// Constraint on the position of a robot point
    ///
    /// The value of the function is defined as the vector from a local
    /// point in the joint frame to a target point in the world frame expressed
    /// in some frame:
    /**
     *  \f{eqnarray*}
     *  \mathbf{e} &=& S B^T \left(\mathbf{x}^* -
     *  M(\mathbf{q})\mathbf{x}\right)
     *  = S B^T \left(\mathbf{x}^* - R(\mathbf{q})\mathbf{x} -
     *    \mathbf{t}(\mathbf{q})\right)\\
     *  \mathbf{\dot{e}} &=& S B^T \left( -
     *  \left[\omega\right]_{\times}R(\mathbf{q})\mathbf{x}
     *                       - \mathbf{\dot{t}}\right)\\
     *                   &=& S B^T \left(
     *  \left[R(\mathbf{q})\mathbf{x}\right]_{\times}\omega
     *                       - \mathbf{\dot{t}}\right)\\
     *                   &=& S B^T \left(
     *  \left[R(\mathbf{q})\mathbf{x}\right]_{\times} J_{joint}^{\omega}
     *                       - J_{joint}^{\mathbf{v}}\right)\mathbf{\dot{q}}
     *  \f}
    **/
    /// where
    ///   \li \f$S\f$ is a selection matrix (diagonal matrix with 1 or 0)
    ///       defining whether each dof is constrained,
    ///   \li \f$B\f$ is a constant rotation matrix representing the basis in
    ///       which the error is expressed
    ///   \li \f$\mathbf{x}^*\f$ is the target point in world frame,
    ///   \li \f$M(\mathbf{q})\f$ is the position of the constrained joint
    ///       as an homogeneous matrix, and
    ///   \li \f$\mathbf{x}\f$ is the position of the local point in the joint
    ///       frame.
    /// \note The constant rotation matrix \f$R\f$ combined with
    ///  selection of degrees of freedom enable users to define
    ///  constraints in a plane or along a line not necessarily
    ///  aligned with the world reference frame axes.
    class HPP_CONSTRAINTS_DLLAPI Position : public DifferentiableFunction
    {
    public:
      /// Identity matrix of size 3
      static matrix3_t I3;
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      /// Return a shared pointer to a new instance
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraint applies to,
      /// \param joint the joint the position of which is constrained,
      /// \param pointInLocalFrame point in local frame,
      /// \param targetInGlobalFrame target point in globalframe,
      /// \param rotation frame in which the error is expressed,
      /// \param mask which component of the error vector to take into
      ///        account.
      static PositionPtr_t create (const std::string& name,
                                   const DevicePtr_t& robot,
                                   const JointPtr_t& joint,
                                   const vector3_t& pointInLocalFrame,
                                   const vector3_t& targetInGlobalFrame,
                                   const matrix3_t& rotation = matrix3_t::getIdentity (),
                                   std::vector <bool> mask = boost::assign::list_of (true)(true)(true));

      /// Return a shared pointer to a new instance
      ///
      /// \param robot the robot the constraint applies to,
      /// \param joint the joint the position of which is constrained,
      /// \param pointInLocalFrame point in local frame,
      /// \param targetInGlobalFrame target point in globalframe,
      /// \param rotation frame in which the error is expressed,
      /// \param mask which component of the error vector to take into
      ///        account.
      static PositionPtr_t create (const DevicePtr_t& robot,
				   const JointPtr_t& joint,
				   const vector3_t& pointInLocalFrame,
				   const vector3_t& targetInGlobalFrame,
				   const matrix3_t& rotation =
				   matrix3_t::getIdentity (),
				   std::vector <bool> mask =
				   boost::assign::list_of (true)(true)(true));

      virtual ~Position () throw () {}

      /// Set position of target point in global frame
      void reference (const vector3_t& reference)
      {
	targetInGlobalFrame_ = reference;
      }

      /// Get position of target point in global frame
      const vector3_t& reference () const
      {
	return targetInGlobalFrame_;
      }
      /// Constructor
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraint applies to,
      /// \param joint the joint the position of which is constrained,
      /// \param pointInLocalFrame point in local frame,
      /// \param targetInGlobalFrame target point in globalframe.
      Position (const std::string& name, const DevicePtr_t& robot,
               const JointPtr_t& joint, const vector3_t& pointInLocalFrame,
		const vector3_t& targetInGlobalFrame, const matrix3_t& rotation,
		std::vector <bool> mask);
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
      vector3_t pointInLocalFrame_;
      vector3_t targetInGlobalFrame_;
      // Case where all degrees of freedom are seleted and rotation is identity
      bool nominalCase_;
      matrix_t SBT_;
      mutable vector3_t p_;
      mutable eigen::matrix3_t cross_;
      mutable vector_t result_;
    }; // class Position
    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_POSITION_HH
