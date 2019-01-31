//
// Copyright (c) 2016 CNRS
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

#ifndef HPP_CONSTRAINTS_GENERIC_TRANSFORMATION_HH
# define HPP_CONSTRAINTS_GENERIC_TRANSFORMATION_HH

# include <pinocchio/spatial/se3.hpp>

# include <hpp/pinocchio/joint.hh>

# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/config.hh>
# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/matrix-view.hh>

namespace hpp {
  namespace constraints {
    /// \cond DEVEL
    template <bool rel> struct GenericTransformationModel
    {
      JointConstPtr_t joint2;
      bool R1isID, R2isID, t1isZero, t2isZero;
      Transform3f F1inJ1, F2inJ2;
      bool fullPos, fullOri;
      size_type rowOri;
      const size_type cols;

      inline JointConstPtr_t getJoint1() const { return JointConstPtr_t(); }
      inline void setJoint1(const JointConstPtr_t&) {}
      GenericTransformationModel (const size_type nCols) :
        joint2(), R1isID(true), R2isID(true), t1isZero(true), t2isZero(true),
        fullPos(false), fullOri(false), cols (nCols)
      { F1inJ1.setIdentity(); F2inJ2.setIdentity(); }
      void checkIsIdentity1() {
        R1isID = F1inJ1.rotation().isIdentity();
        t1isZero = F1inJ1.translation().isZero();
      }
      void checkIsIdentity2() {
        R2isID = F2inJ2.rotation().isIdentity();
        t2isZero = F2inJ2.translation().isZero();
      }
    };
    template <> struct GenericTransformationModel<true> :
      GenericTransformationModel<false>
    {
      JointConstPtr_t joint1;
      inline JointConstPtr_t getJoint1() const { return joint1; }
      void setJoint1(const JointConstPtr_t& j);
      GenericTransformationModel (const size_type nCols) :
        GenericTransformationModel<false>(nCols), joint1() {}
    };
    /// \endcond DEVEL

    /// \addtogroup constraints
    /// \{

    /** GenericTransformation class encapsulates 6 possible differentiable
     *  functions: Position, Orientation, Transformation and their relative
     *  versions RelativePosition, RelativeOrientation, RelativeTransformation.
     *
     *  These functions compute the position of frame
     *  GenericTransformation::frame2InJoint2 in joint
     *  GenericTransformation::joint2 frame, in the frame
     *  GenericTransformation::frame1InJoint1 in GenericTransformation::joint1
     *  frame. For absolute functions, GenericTransformation::joint1 is
     *  NULL and joint1 frame is the world frame.
     *
     *  The value of the RelativeTransformation function is a 6-dimensional
     *  vector. The 3 first coordinates are the position of the center of the
     *  second frame expressed in the first frame.
     *  The 3 last coordinates are the log of the orientation of frame 2 in
     *  frame 1.
     *
     *  \f{equation*}
     *  f (\mathbf{q}) = \left(\begin{array}{c}
     *  \mathbf{translation}\left(T_{1/J_1}^T T_1^T T_2 T_{2/J_2}\right)\\
     *  \log \left((R_1 R_{1/J_1})^T R_2 R_{2/J_2}\right) \end{array}\right)
     *  \f}
     *
     *  The Jacobian is given by
     *
     *  \f{equation*}
     *  \left(\begin{array}{c}
     *  (R_1 R_{1/J_1})^T \left(\left[R_2 t_{2/J_2} + t_2 - t_1\right]_{\times}
     *  R_1 J_{1\,\omega} - \left[R_2 t_{2/J_2}\right]_{\times} R_2 J_{2\,\omega} +
     *  R_2 J_{2\,\mathbf{v}} - R_1 J_{1\,\mathbf{v}}\right) \\
     *  J_{log}\left((R_1 R_{1/J_1})^T R_2 R_{2/J_2}\right)(R_1 R_{1/J_1})^T
     *  (R_2 J_{2\,\omega} - R_1 J_{1\,\omega})
     *  \end{array}\right)
     *  \f}
    */
    template <int _Options>
    class HPP_CONSTRAINTS_DLLAPI GenericTransformation :
      public DifferentiableFunction
    {
    public:
      typedef boost::shared_ptr <GenericTransformation> Ptr_t;
      typedef boost::weak_ptr <GenericTransformation> WkPtr_t;

      enum {
        IsRelative         = _Options & RelativeBit,
        ComputeOrientation = _Options & OrientationBit,
        ComputePosition    = _Options & PositionBit,
        OutputSE3          = _Options & OutputSE3Bit,
        IsPosition         = ComputePosition  && !ComputeOrientation,
        IsOrientation      = !ComputePosition && ComputeOrientation,
        IsTransform        = ComputePosition  && ComputeOrientation,
        ValueSize          = (ComputePosition?3:0) + (ComputeOrientation?(OutputSE3?4:3):0),
        DerSize            = (ComputePosition?3:0) + (ComputeOrientation ?3:0)
      };

      /// \cond
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      /// \endcond

      /// Object builder for absolute functions.
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      /// \param joint2 the second joint the transformation of which is
      ///               constrained,
      /// \param reference desired relative transformation
      ///        \f$T_2(\mathbf{q})\f$ between the joints.
      /// \param mask which component of the error vector to take into
      ///        account.
      static Ptr_t create (const std::string& name, const DevicePtr_t& robot,
          const JointConstPtr_t& joint2, const Transform3f& reference,
          std::vector <bool> mask = std::vector<bool>(DerSize,true));

      /// Object builder for absolute functions.
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      /// \param joint2 the second joint the transformation of which is
      ///               constrained,
      /// \param frame2 position of a fixed frame in joint 2,
      /// \param frame1 position of a fixed frame in the world,
      /// \param mask vector of 6 boolean defining which coordinates of the
      ///        error vector to take into account.
      ///
      /// \note For Position, the rotation part of frame1 defines the
      ///       frame in which the error is expressed and the rotation of frame2
      ///       has no effect.
      static Ptr_t create (const std::string& name, const DevicePtr_t& robot,
          /* World frame          */ const JointConstPtr_t& joint2,
          const Transform3f& frame2, const Transform3f& frame1,
         std::vector <bool> mask = std::vector<bool>(DerSize,true));

      /// Object builder for relative functions.
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
      static Ptr_t create (const std::string& name, const DevicePtr_t& robot,
          const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
          const Transform3f& reference,
          std::vector <bool> mask = std::vector<bool>(DerSize,true));

      /// Object builder for relative functions.
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
      ///
      /// \note For RelativePosition, the rotation part of frame1 defines the
      ///       frame in which the error is expressed and the rotation of frame2
      ///       has no effect.
      static Ptr_t create (const std::string& name, const DevicePtr_t& robot,
	 const JointConstPtr_t& joint1,  const JointConstPtr_t& joint2,
	 const Transform3f& frame1, const Transform3f& frame2,
         std::vector <bool> mask = std::vector<bool>(DerSize,true));

      virtual ~GenericTransformation () {}

      /// Set desired relative transformation of joint2 in joint1
      ///
      inline void reference (const Transform3f& reference)
      {
	m_.F1inJ1 = reference;
        m_.checkIsIdentity1();
	m_.F2inJ2.setIdentity ();
        m_.checkIsIdentity2();
      }

      /// Get desired relative orientation
      inline Transform3f reference () const
      {
	return m_.F1inJ1.actInv(m_.F2inJ2);
      }

      /// Set joint 1
      inline void joint1 (const JointConstPtr_t& joint) {
        // static_assert(IsRelative);
	m_.setJoint1(joint);
        computeActiveParams();
	assert (!joint || joint->robot () == robot_);
      }

      /// Get joint 1
      inline JointConstPtr_t joint1 () const {
	return m_.getJoint1();
      }

      /// Set joint 2
      inline void joint2 (const JointConstPtr_t& joint) {
	m_.joint2 = joint;
        computeActiveParams();
	assert (!joint || (joint->index() > 0 && joint->robot () == robot_));
      }

      /// Get joint 2
      inline JointConstPtr_t joint2 () const {
	return m_.joint2;
      }

      /// Set position of frame 1 in joint 1
      inline void frame1InJoint1 (const Transform3f& M) {
	m_.F1inJ1 = M;
        m_.checkIsIdentity1();
      }
      /// Get position of frame 1 in joint 1
      inline const Transform3f& frame1InJoint1 () const {
	return m_.F1inJ1;
      }
      /// Set position of frame 2 in joint 2
      inline void frame2InJoint2 (const Transform3f& M) {
	m_.F2inJ2 = M;
        m_.checkIsIdentity2();
      }
      /// Get position of frame 2 in joint 2
      inline const Transform3f& frame2InJoint2 () const {
	return m_.F2inJ2;
      }

      virtual std::ostream& print (std::ostream& o) const;

      ///Constructor
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      ///        \f$T_1(\mathbf{q})^{-1} T_2(\mathbf{q})\f$ between the joints.
      /// \param mask vector of 6 boolean defining which coordinates of the
      ///        error vector to take into account.
      GenericTransformation (const std::string& name,
                              const DevicePtr_t& robot,
                              std::vector <bool> mask);

    protected:
      void init (const WkPtr_t& self)
      {
        self_ = self;
        computeActiveParams();
      }

      /// Compute value of error
      ///
      /// \param argument configuration of the robot,
      /// \retval result error vector
      virtual void impl_compute	(LiegroupElementRef result,
				 ConfigurationIn_t argument) const throw ();
      virtual void impl_jacobian (matrixOut_t jacobian,
				  ConfigurationIn_t arg) const throw ();
    private:
      void computeActiveParams ();
      DevicePtr_t robot_;
      GenericTransformationModel<IsRelative> m_;
      Eigen::RowBlockIndices Vindices_;
      const std::vector <bool> mask_;
      WkPtr_t self_;
    }; // class GenericTransformation
    /// \}
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_GENERIC_TRANSFORMATION_HH
