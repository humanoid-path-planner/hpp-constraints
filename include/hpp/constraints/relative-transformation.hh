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

#ifndef HPP_RELATIVE_TRANSFORMATION_HH
# define HPP_RELATIVE_TRANSFORMATION_HH

# include <boost/assign/list_of.hpp>
# include <hpp/fcl/math/transform.h>
# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/relative-orientation.hh>
# include <hpp/constraints/relative-position.hh>
# include <hpp/constraints/config.hh>
# include <hpp/constraints/fwd.hh>

namespace hpp {
  namespace constraints {
    /// \addtogroup constraints
    /// \{

    /** Relative transformation of two fixed frames in robot joints

        The value of the function is a 6-dimensional vector. The 3 first
        coordinates are the position of the center of the second frame expressed
        in the first frame. The 3 last coordinates are the log of the
        orientation of frame 2 in frame 1.

        \f{equation*}
	f (\mathbf{q}) = \left(\begin{array}{c}
	\mathbf{translation}\left(T_{1/J_1}^T T_1^T T_2 T_{2/J_2}\right)\\
	\log ((R_1 R_{1/J_1})^T R_2 R_{2/J_2}) \end{array}\right)
        \f}

	The Jacobian is given by

	\f{equation*}
	\left(\begin{array}{c}
	(R_1 R_{1/J_1})^T (\left[R_2 t_{2/J_2} + t_2 - t_1\right]_{\times}
	J_{1\,\omega} - \left[R_2 t_{2/J_2}\right]_{\times} J_{2\,\omega} +
	J_{2\,\mathbf{v}} - J_{1\,\mathbf{v}}) \\
	J_{log}((R_1 R_{1/J_1})^T R_2 R_{2/J_2})(R_1 R_{1/J_1})^T
	(J_{2\,\omega} - J_{1\,\omega})
	\end{array}\right)
	\f}
    */
    class HPP_CONSTRAINTS_DLLAPI RelativeTransformation :
      public DifferentiableFunction
    {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

      /// Return a shared pointer to a new instance
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      /// \param joint1 the first joint the transformation of which is
      ///               constrained,
      /// \param joint2 the second joint the transformation of which is
      ///               constrained,
      /// \param reference desired relative transformation
      ///        \f$T_1(\mathbf{q})^{-1} T_2(\mathbf{q})\f$ between the joints.
      /// \param mask which component of the error vector to take into
      ///        account.
      static RelativeTransformationPtr_t create (
          const std::string& name,
          const DevicePtr_t& robot,
          const JointPtr_t& joint1,
          const JointPtr_t& joint2,
          const Transform3f& reference,
          std::vector <bool> mask =
            boost::assign::list_of (true)(true)(true)
                                   (true)(true)(true));

      /// Return a shared pointer to a new instance
      ///
      /// \param robot the robot the constraints is applied to,
      /// \param joint1 the first joint the transformation of which is
      ///               constrained,
      /// \param joint2 the second joint the transformation of which is
      ///               constrained,
      /// \param reference desired relative transformation
      ///        \f$T_1(\mathbf{q})^{-1} T_2(\mathbf{q})\f$ between the joints.
      /// \param mask which component of the error vector to take into
      ///        account.
      static RelativeTransformationPtr_t create
	(const DevicePtr_t& robot, const JointPtr_t& joint1,
	 const JointPtr_t& joint2, const Transform3f& reference,
	 std::vector <bool> mask = boost::assign::list_of (true)(true)(true)
	 (true)(true)(true));

      /// Return a shared pointer to a new instance
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      /// \param joint1 the first joint the transformation of which is
      ///               constrained,
      /// \param joint2 the second joint the transformation of which is
      ///               constrained,
      /// \param frame1 position of a fixed frame in joint 1,
      /// \param frame2 position of a fixed frame in joint 2,
      /// \param mask vector of 6 boolean defining which coordinates of the
      ///        error vector to take into account.
      /// \note if joint1 is 0x0, joint 1 frame is considered to be the global
      ///       frame.
      static RelativeTransformationPtr_t create
	(const std::string& name, const DevicePtr_t&,
	 const JointPtr_t& joint1, const JointPtr_t& joint2,
	 const Transform3f& frame1, const Transform3f& frame2,
	 std::vector <bool> mask = boost::assign::list_of (true)(true)(true)
	 (true)(true)(true));

      virtual ~RelativeTransformation () throw () {}
      /// Set desired relative transformation of joint2 in joint1
      ///
      void reference (const Transform3f& reference)
      {
	F1inJ1_ = reference;
	F2inJ2_.setIdentity ();
      }

      /// Get desired relative orientation
      Transform3f reference () const
      {
	return F1inJ1_ * inverse (F2inJ2_);
      }

      /// Set joint 1
      void joint1 (const JointPtr_t& joint) {
	joint1_ = joint;
	assert (!joint || joint->robot () == robot_);
      }

      /// Get joint 1
      JointPtr_t joint1 () {
	return joint1_;
      }

      /// Set joint 2
      void joint2 (const JointPtr_t& joint) {
	joint2_ = joint;
	assert (!joint || joint->robot () == robot_);
      }

      /// Get joint 2
      JointPtr_t joint2 () {
	return joint2_;
      }

      /// Set position of frame 1 in joint 1
      void frame1inJoint1 (const Transform3f& M) {
	F1inJ1_ = M;
      }
      /// Get position of frame 1 in joint 1
      const Transform3f& frame1inJoint1 () const {
	return F1inJ1_;
      }
      /// Set position of frame 2 in joint 2
      void frame2inJoint2 (const Transform3f& M) {
	F2inJ2_ = M;
      }
      /// Get position of frame 2 in joint 2
      const Transform3f& frame2inJoint2 () const {
	return F2inJ2_;
      }
    protected:
      ///Constructor
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      /// \param joint1 the first joint the transformation of which is
      ///               constrained,
      /// \param joint2 the second joint the transformation of which is
      ///               constrained,
      /// \param reference desired relative transformation
      ///        \f$T_1(\mathbf{q})^{-1} T_2(\mathbf{q})\f$ between the joints.
      /// \param mask vector of 6 boolean defining which coordinates of the
      ///        error vector to take into account.
      RelativeTransformation (const std::string& name,
                              const DevicePtr_t&,
                              const JointPtr_t& joint1,
                              const JointPtr_t& joint2,
                              const Transform3f& reference,
                              std::vector <bool> mask);

      ///Constructor
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      /// \param joint1 the first joint the transformation of which is
      ///               constrained,
      /// \param joint2 the second joint the transformation of which is
      ///               constrained,
      /// \param frame1 position of a fixed frame in joint 1,
      /// \param frame2 position of a fixed frame in joint 2,
      /// \param mask vector of 6 boolean defining which coordinates of the
      ///        error vector to take into account.
      RelativeTransformation (const std::string& name, const DevicePtr_t&,
                              const JointPtr_t& joint1,
			      const JointPtr_t& joint2,
                              const Transform3f& frame1,
                              const Transform3f& frame2,
                              std::vector <bool> mask);
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
      void computeError (ConfigurationIn_t argument) const;
      DevicePtr_t robot_;
      JointPtr_t joint1_;
      JointPtr_t joint2_;
      Transform3f F1inJ1_;
      Transform3f F2inJ2_;
      const std::vector <bool> mask_;
      mutable value_type theta_;
      mutable vector_t value_;
      mutable matrix_t jacobian_;
      mutable eigen::matrix3_t cross1_, cross2_, Jlog_;
      mutable Configuration_t latestArgument_;
    }; // class RelativeTransformation
    /// \}
  } // namespace constraints
} // namespace hpp
#endif // HPP_RELATIVE_TRANSFORMATION_HH
