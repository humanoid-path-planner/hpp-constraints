//
// Copyright (c) 2016 CNRS
// Authors: Joseph Mirabel
//
//

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#ifndef HPP_CONSTRAINTS_GENERIC_TRANSFORMATION_HH
# define HPP_CONSTRAINTS_GENERIC_TRANSFORMATION_HH

# include <pinocchio/spatial/se3.hpp>

# include <hpp/util/serialization-fwd.hh>

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

    /** GenericTransformation class encapsulates 10 possible differentiable
     *  functions: Position, Orientation, OrientationSO3, Transformation,
     *  TransformationR3xSO3 and their relative
     *  versions RelativeTransformationR3xSO3, RelativePosition,
     *  RelativeOrientation, RelativeOrientationSO3, RelativeTransformation.
     *
     *  These functions compute the position of frame
     *  GenericTransformation::frame2InJoint2 in joint
     *  GenericTransformation::joint2 frame, in the frame
     *  GenericTransformation::frame1InJoint1 in GenericTransformation::joint1
     *  frame. For absolute functions, GenericTransformation::joint1 is
     *  NULL and joint1 frame is the world frame.
     *
     *  The value of the RelativeTransformation function is
     *  a 6-dimensional
     *  vector. The 3 first coordinates are the position of the center of the
     *  second frame expressed in the first frame.
     *  The 3 last coordinates are the log of the orientation of frame 2 in
     *  frame 1.
     *  The values of RelativePosition and RelativeOrientation are respectively
     *  the 3 first and the 3 last components of the above 6 dimensional vector.
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
     *
     *  For RelativeTransformationR3xSO3, OrientationSO3, the 3 last components
     *  are replaced by a quaternion.
    */
    template <int _Options>
    class HPP_CONSTRAINTS_DLLAPI GenericTransformation :
      public DifferentiableFunction
    {
    public:
      typedef shared_ptr <GenericTransformation> Ptr_t;
      typedef weak_ptr <GenericTransformation> WkPtr_t;
#if __cplusplus >= 201103L
      static constexpr bool
        IsRelative         = _Options & RelativeBit,
        ComputeOrientation = _Options & OrientationBit,
        ComputePosition    = _Options & PositionBit,
        OutputR3xSO3          = _Options & OutputR3xSO3Bit,
        IsPosition         = ComputePosition  && !ComputeOrientation,
        IsOrientation      = !ComputePosition && ComputeOrientation,
        IsTransform        = ComputePosition  && ComputeOrientation;
      static constexpr int
        ValueSize          = (ComputePosition?3:0) +
	(ComputeOrientation?(OutputR3xSO3?4:3):0),
        DerSize            = (ComputePosition?3:0) + (ComputeOrientation ?3:0);
#else
      enum {
        IsRelative         = _Options & RelativeBit,
        ComputeOrientation = _Options & OrientationBit,
        ComputePosition    = _Options & PositionBit,
        OutputR3xSO3          = _Options & OutputR3xSO3Bit,
        IsPosition         = ComputePosition  && !ComputeOrientation,
        IsOrientation      = !ComputePosition && ComputeOrientation,
        IsTransform        = ComputePosition  && ComputeOrientation,
        ValueSize          = (ComputePosition?3:0) +
	(ComputeOrientation?(OutputR3xSO3?4:3):0),
        DerSize            = (ComputePosition?3:0) + (ComputeOrientation ?3:0)
        };
#endif

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
				 ConfigurationIn_t argument) const;
      virtual void impl_jacobian (matrixOut_t jacobian,
				  ConfigurationIn_t arg) const;

      bool isEqual(const DifferentiableFunction& other) const {
        const GenericTransformation& castother = dynamic_cast<const GenericTransformation&>(other);
        if (!DifferentiableFunction::isEqual(other))
          return false;
        if (robot_ != castother.robot_)
          return false;
        if (mask_ != castother.mask_)
          return false;
        
        return true;
      }

    private:
      void computeActiveParams ();
      DevicePtr_t robot_;
      GenericTransformationModel<IsRelative> m_;
      Eigen::RowBlockIndices Vindices_;
      const std::vector <bool> mask_;
      WkPtr_t self_;

      GenericTransformation() : m_ (0) {}
      HPP_SERIALIZABLE();
    }; // class GenericTransformation
    /// \}
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_GENERIC_TRANSFORMATION_HH
