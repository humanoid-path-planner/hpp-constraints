///
/// Copyright (c) 2014 CNRS
/// Authors: Florent Lamiraux
///
///
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

#ifndef HPP_CONSTRAINTS_ORIENTATION_HH
# define HPP_CONSTRAINTS_ORIENTATION_HH

# include <boost/assign/list_of.hpp>
# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/config.hh>
# include <hpp/constraints/fwd.hh>

namespace hpp {
  namespace constraints {
    /// \addtogroup constraints
    /// \{

    /// Constraint on the orientation of a robot joint
    ///
    /// The value of the function is defined as the rotation vector of the
    /// reference orientation expressed in the frame of the current
    /// orientation.
    /**
     *  \f{eqnarray*}
     *  \mathbf{r} (\mathbf{q}) &=& \log \left(R(\mathbf{q})^TR^*)\right)\\
     *  \mathbf{\dot{r}} &=& -J_{\log}(\left(R(\mathbf{q})^TR^*)\right))
     *  R(\mathbf{q})^T J^{\omega}_{joint}(\mathbf{q}) \mathbf{\dot{q}}\\
     *  \f}
    **/
    /// where
    /**
     *  \f{eqnarray*}
     *  J_{\log} (R) &=& \frac{\|\mathbf{r}\|\sin\|\mathbf{r}\|}
     *                   {2(1-\cos\|\mathbf{r}\|)}I_3
     *                   -\frac{1}{2}[\mathbf{r}]_{\times} +
     *                   \left(\frac{1}{\|\mathbf{r}\|^2} -
     *                   \frac{\sin\|\mathbf{r}\|}{2\|\mathbf{r}\|
     *                   (1-\cos\|\mathbf{r}\|)}\right)
     *                   \mathbf{r}\mathbf{r}^T
     *  \f}
    **/
    class HPP_CONSTRAINTS_DLLAPI Orientation : public DifferentiableFunction
    {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      /// Return a shared pointer to a new instance
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      /// \param joint the joint the orientation of which is constrained
      /// \param reference reference orientation of the joint,
      /// \param mask which component of the error vector to take into
      ///        account.
      static OrientationPtr_t create (const std::string& name,
                                      const DevicePtr_t& robot,
                                      const JointPtr_t& joint,
                                      const matrix3_t& reference,
                                      std::vector <bool> mask =
                                      boost::assign::list_of (true)(true)(true));

      /// Return a shared pointer to a new instance
      ///
      /// \param robot the robot the constraints is applied to,
      /// \param joint the joint the orientation of which is constrained
      /// \param reference reference orientation of the joint,
      /// \param mask which component of the error vector to take into
      ///        account.
      static OrientationPtr_t create (const DevicePtr_t& robot,
				      const JointPtr_t& joint,
				      const matrix3_t& reference,
				      std::vector <bool> mask =
				      boost::assign::list_of (true)(true)(true));
      virtual ~Orientation () throw () {}
      /// Set desired orientation as a rotation matrix
      void reference (const matrix3_t& reference)
      {
	reference_ = reference;
      }
      /// Get desired orientation
      const matrix3_t& reference () const
      {
	return reference_;
      }

      ///Constructor
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      /// \param joint the joint the orientation of which is constrained
      /// \param reference reference orientation of the joint,
      /// \param mask which component of the error vector to take into
      ///        account.
      Orientation (const std::string& name, const DevicePtr_t&,
                   const JointPtr_t& joint, const matrix3_t& reference,
                   std::vector <bool> mask);

      /// Get size of error with respect to mask.
      static size_type size (std::vector<bool> mask);
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
      void computeError (vectorOut_t result, ConfigurationIn_t argument,
			 double& theta, std::vector<bool> mask) const;
      DevicePtr_t robot_;
      JointPtr_t joint_;
      matrix3_t reference_;
      matrix3_t ref_;
      std::vector <bool> mask_;
      mutable matrix3_t Rerror_;
      mutable vector_t r_;
      mutable eigen::matrix3_t Jlog_;
      mutable matrix_t jacobian_;
    }; // class Orientation
    /// \}
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_ORIENTATION_HH
