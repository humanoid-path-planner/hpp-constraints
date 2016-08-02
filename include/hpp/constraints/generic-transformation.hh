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

namespace hpp {
  namespace constraints {
    /// \cond DEVEL
    template <bool rel> struct GenericTransformationJointData
    {
      JointPtr_t joint2;
      bool R1isID, R2isID, t1isZero, t2isZero;
      Transform3f F1inJ1, F2inJ2;
      inline JointPtr_t getJoint1() const { return JointPtr_t(); }
      inline void setJoint1(const JointPtr_t&) {}
      const JointJacobian_t& J2 () const { return joint2->jacobian(); }
      const Transform3f& M2 () const { return joint2->currentTransformation(); }
      const vector3_t& t2 () const { return joint2->currentTransformation().translation(); }
      const matrix3_t& R2 () const { return joint2->currentTransformation().rotation(); }
      GenericTransformationJointData () :
        joint2(), R1isID(true), R2isID(true), t1isZero(true), t2isZero(true)
      { F1inJ1.setIdentity(); F2inJ2.setIdentity(); }
    };
    template <> struct GenericTransformationJointData<true> :
      GenericTransformationJointData<false>
    {
      JointPtr_t joint1;
      inline JointPtr_t getJoint1() const { return joint1; }
      inline void setJoint1(const JointPtr_t& j) { joint1 = j; }
      const JointJacobian_t& J1 () const { return joint1->jacobian(); }
      const Transform3f& M1 () const { return joint1->currentTransformation(); }
      const matrix3_t& R1 () const { return joint1->currentTransformation().rotation(); }
      const vector3_t& t1 () const { return joint1->currentTransformation().translation(); }
      GenericTransformationJointData () :
        GenericTransformationJointData<false>(), joint1() {}
    };
    template <bool rel> struct GenericTransformationOriData {};
    template <> struct GenericTransformationOriData<true>
    {
      mutable Transform3f M;
      mutable value_type theta;
      mutable eigen::matrix3_t JlogXTR1inJ1;
    };
    /// This class contains the data of the GenericTransformation class.
    template <bool rel, bool pos, bool ori> struct GenericTransformationData :
      GenericTransformationJointData<rel>,
      GenericTransformationOriData<ori>
    {
      enum { NbRows = (pos?3:0)+(ori?3:0) };
      typedef Eigen::Matrix<value_type, NbRows, 1> ValueType;
      typedef Eigen::Matrix<value_type, NbRows, Eigen::Dynamic> JacobianType;
      bool fullPos, fullOri;
      size_type rowOri;
      const size_type cols;
      mutable ValueType value;
      mutable JacobianType jacobian;
      mutable Eigen::Matrix<value_type, 3, Eigen::Dynamic> tmpJac;
      mutable eigen::vector3_t cross1, cross2;
      GenericTransformationData (const size_type nCols) :
        GenericTransformationJointData<rel>(),
        GenericTransformationOriData<ori> (),
        fullPos(false), fullOri(false), cols (nCols),
        jacobian((int)NbRows, cols)
      { cross1.setZero(); cross2.setZero(); }
      void checkIsIdentity1() {
        this->R1isID = this->F1inJ1.rotation().isIdentity(); this->t1isZero = this->F1inJ1.translation().isZero();
      }
      void checkIsIdentity2() {
        this->R2isID = this->F2inJ2.rotation().isIdentity(); this->t2isZero = this->F2inJ2.translation().isZero();
        if (this->t2isZero) cross2.setZero();
      }
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
     *  \log ((R_1 R_{1/J_1})^T R_2 R_{2/J_2}) \end{array}\right)
     *  \f}
     *
     *  The Jacobian is given by
     *
     *  \f{equation*}
     *  \left(\begin{array}{c}
     *  (R_1 R_{1/J_1})^T (\left[R_2 t_{2/J_2} + t_2 - t_1\right]_{\times}
     *  J_{1\,\omega} - \left[R_2 t_{2/J_2}\right]_{\times} J_{2\,\omega} +
     *  J_{2\,\mathbf{v}} - J_{1\,\mathbf{v}}) \\
     *  J_{log}((R_1 R_{1/J_1})^T R_2 R_{2/J_2})(R_1 R_{1/J_1})^T
     *  (J_{2\,\omega} - J_{1\,\omega})
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
        IsPosition         = ComputePosition  && !ComputeOrientation,
        IsOrientation      = !ComputePosition && ComputeOrientation,
        IsTransform        = ComputePosition  && ComputeOrientation,
        ValueSize          = (ComputePosition?3:0) + (ComputeOrientation?3:0)
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
          const JointPtr_t& joint2, const Transform3f& reference,
          std::vector <bool> mask = std::vector<bool>(ValueSize,true));

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
          /* World frame          */ const JointPtr_t& joint2,
          const Transform3f& frame2, const Transform3f& frame1,
         std::vector <bool> mask = std::vector<bool>(ValueSize,true));

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
          const JointPtr_t& joint1, const JointPtr_t& joint2,
          const Transform3f& reference,
          std::vector <bool> mask = std::vector<bool>(ValueSize,true));

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
	 const JointPtr_t& joint1,  const JointPtr_t& joint2,
	 const Transform3f& frame1, const Transform3f& frame2,
         std::vector <bool> mask = std::vector<bool>(ValueSize,true));

      virtual ~GenericTransformation () {}

      /// Set desired relative transformation of joint2 in joint1
      ///
      inline void reference (const Transform3f& reference)
      {
	d_.F1inJ1 = reference;
        d_.checkIsIdentity1();
	d_.F2inJ2.setIdentity ();
        d_.checkIsIdentity2();
      }

      /// Get desired relative orientation
      inline Transform3f reference () const
      {
	return d_.F1inJ1.actInv(d_.F2inJ2);
      }

      /// Set joint 1
      inline void joint1 (const JointPtr_t& joint) {
        // static_assert(IsRelative);
	d_.setJoint1(joint);
	assert (!joint || joint->robot () == robot_);
      }

      /// Get joint 1
      inline JointPtr_t joint1 () {
	return d_.getJoint1();
      }

      /// Set joint 2
      inline void joint2 (const JointPtr_t& joint) {
	d_.joint2 = joint;
	assert (!joint || joint->robot () == robot_);
      }

      /// Get joint 2
      inline JointPtr_t joint2 () {
	return d_.joint2;
      }

      /// Set position of frame 1 in joint 1
      inline void frame1InJoint1 (const Transform3f& M) {
	d_.F1inJ1 = M;
        d_.checkIsIdentity1();
      }
      /// Get position of frame 1 in joint 1
      inline const Transform3f& frame1InJoint1 () const {
	return d_.F1inJ1;
      }
      /// Set position of frame 2 in joint 2
      inline void frame2InJoint2 (const Transform3f& M) {
	d_.F2inJ2 = M;
        d_.checkIsIdentity2();
      }
      /// Get position of frame 2 in joint 2
      inline const Transform3f& frame2InJoint2 () const {
	return d_.F2inJ2;
      }

      ///Constructor
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      ///        \f$T_1(\mathbf{q})^{-1} T_2(\mathbf{q})\f$ between the joints.
      /// \param mask vector of 6 boolean defining which coordinates of the
      ///        error vector to take into account.
      GenericTransformation (const std::string& name,
                              const DevicePtr_t&,
                              std::vector <bool> mask);

    protected:
      void init (const WkPtr_t& self)
      {
        self_ = self;
      }

      /// Compute value of error
      ///
      /// \param argument configuration of the robot,
      /// \retval result error vector
      virtual void impl_compute	(vectorOut_t result,
				 ConfigurationIn_t argument) const throw ();
      virtual void impl_jacobian (matrixOut_t jacobian,
				  ConfigurationIn_t arg) const throw ();
    private:
      void computeError (const ConfigurationIn_t& argument) const;
      DevicePtr_t robot_;
      GenericTransformationData<IsRelative,ComputePosition,ComputeOrientation>
        d_;
      const std::vector <bool> mask_;
      WkPtr_t self_;
      mutable Configuration_t latestArgument_;
    }; // class GenericTransformation
    /// \}
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_GENERIC_TRANSFORMATION_HH
