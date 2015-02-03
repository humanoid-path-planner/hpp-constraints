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

    /// Constraint on the relative transformation of two robot joints
    ///
    /// The 3 first coordinates of the error is computed by a RelativePosition
    /// instance where the origin of Joint 2 should coincide with the
    /// translation part of the reference transformation passed to this.
    ///
    /// The 3 last coordinates of the error are computed by a
    /// RelativeOrientation instance.
    ///
    class HPP_CONSTRAINTS_DLLAPI RelativeTransformation :
      public DifferentiableFunction
    {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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
      virtual ~RelativeTransformation () throw () {}
      /// Set desired relative transformation
      void reference (const Transform3f& reference)
      {
	reference_ = reference;
	relativeOrientation_.reference (reference_.getRotation ());
	relativePosition_.pointInJoint1 (reference_.getTranslation ());
      }
      /// Get desired relative orientation
      const Transform3f& reference () const
      {
	return reference_;
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
      RelativeOrientation relativeOrientation_;
      RelativePosition relativePosition_;
      Transform3f reference_;
      /// Size of the translation constraint
      size_type sizeTranslation_;
      /// Size of the orientation constraint
      size_type sizeOrientation_;
    }; // class RelativeTransformation
    /// \}
  } // namespace constraints
} // namespace hpp
#endif // HPP_RELATIVE_TRANSFORMATION_HH
