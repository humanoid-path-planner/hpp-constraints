//
// Copyright (c) 2015 CNRS
// Authors: Florent Lamiraux, Mylene Campana
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

#ifndef HPP_TRANSFORMATION_HH
# define HPP_TRANSFORMATION_HH

# include <boost/assign/list_of.hpp>
# include <hpp/fcl/math/transform.h>
# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/orientation.hh>
# include <hpp/constraints/position.hh>
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

    /// Constraint on the transformation of one robot joint
    ///
    /// The 3 first coordinates of the error is computed by a Position
    /// instance where the origin should coincide with the
    /// translation part of the reference transformation passed to this.
    ///
    /// The 3 last coordinates of the error are computed by an
    /// Orientation instance.
    ///
    class HPP_CONSTRAINTS_DLLAPI Transformation :
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
      /// \param reference desired transformation
      ///        \f$T_1(\mathbf{q})^{-1} T_2(\mathbf{q})\f$ between the frames.
      /// \param mask which component of the error vector to take into
      ///        account.
      static TransformationPtr_t create (
          const std::string& name,
          const DevicePtr_t& robot,
          const JointPtr_t& joint1,
          const Transform3f& reference,
          std::vector <bool> mask =
            boost::assign::list_of (true)(true)(true)
                                   (true)(true)(true));

      /// Return a shared pointer to a new instance
      ///
      /// \param robot the robot the constraints is applied to,
      /// \param joint1 the first joint the transformation of which is
      ///               constrained,
      /// \param reference desired transformation
      ///        \f$T_1(\mathbf{q})^{-1} T_2(\mathbf{q})\f$ between the frames.
      /// \param mask which component of the error vector to take into
      ///        account.
      static TransformationPtr_t create
	(const DevicePtr_t& robot, const JointPtr_t& joint1,
	 const Transform3f& reference,
	 std::vector <bool> mask = boost::assign::list_of (true)(true)(true)
	 (true)(true)(true));
      virtual ~Transformation () throw () {}
      /// Set desired transformation
      void reference (const Transform3f& reference)
      {
	reference_ = reference;
	orientation_.reference (reference_.getRotation ());
	// set position TargetInGobalFrame
	position_.reference (reference_.getTranslation ());
      }
      /// Get desired orientation
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
      /// \param reference desired transformation
      ///        \f$T_1(\mathbf{q})^{-1} T_2(\mathbf{q})\f$ between the frames.
      /// \param mask vector of 6 boolean defining which coordinates of the
      ///        error vector to take into account.
      Transformation (const std::string& name,
                              const DevicePtr_t&,
                              const JointPtr_t& joint1,
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
      Orientation orientation_;
      Position position_;
      Transform3f reference_;
      /// Size of the translation constraint
      size_type sizeTranslation_;
      /// Size of the orientation constraint
      size_type sizeOrientation_;
    }; // class Transformation
    /// \}
  } // namespace constraints
} // namespace hpp
#endif // HPP_TRANSFORMATION_HH
