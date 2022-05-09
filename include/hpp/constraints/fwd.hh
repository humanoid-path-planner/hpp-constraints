///
/// Copyright (c) 2014 CNRS
/// Authors: Florent Lamiraux
///
///

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

#ifndef HPP_CONSTRAINTS_FWD_HH
# define HPP_CONSTRAINTS_FWD_HH

# include <list>
# include <hpp/pinocchio/fwd.hh>
# include <hpp/constraints/deprecated.hh>

namespace Eigen {
  struct BlockIndex;
} // namespace Eigen

namespace hpp {
  namespace constraints {
    HPP_PREDEF_CLASS (DifferentiableFunction);
    HPP_PREDEF_CLASS (DifferentiableFunctionSet);
    HPP_PREDEF_CLASS (ActiveSetDifferentiableFunction);
    typedef pinocchio::size_type size_type;
    typedef pinocchio::value_type value_type;
    typedef pinocchio::JointPtr_t JointPtr_t;
    typedef pinocchio::JointConstPtr_t JointConstPtr_t;
    typedef pinocchio::Joint Joint;
    typedef pinocchio::vector3_t vector3_t;
    typedef pinocchio::matrix3_t matrix3_t;
    typedef Eigen::Matrix<value_type, 6, 6> matrix6_t;
    typedef Eigen::Matrix<value_type, 8, 1> vector8_t;
    typedef pinocchio::matrix_t matrix_t;
    typedef Eigen::Ref <const matrix_t> matrixIn_t;
    typedef Eigen::Ref <matrix_t> matrixOut_t;
    typedef pinocchio::vector_t vector_t;
    typedef pinocchio::vectorIn_t vectorIn_t;
    typedef pinocchio::vectorOut_t vectorOut_t;
    typedef pinocchio::ComJacobian_t ComJacobian_t;
    typedef pinocchio::JointJacobian_t JointJacobian_t;
    typedef pinocchio::Transform3f Transform3f;
    typedef pinocchio::LiegroupElement LiegroupElement;
    typedef pinocchio::LiegroupElementRef LiegroupElementRef;
    typedef pinocchio::LiegroupElementConstRef LiegroupElementConstRef;
    typedef pinocchio::LiegroupSpace LiegroupSpace;
    typedef pinocchio::LiegroupSpacePtr_t LiegroupSpacePtr_t;
    typedef pinocchio::LiegroupSpaceConstPtr_t LiegroupSpaceConstPtr_t;
    namespace eigen {
      typedef Eigen::Matrix <value_type, 3, 3> matrix3_t;
      typedef Eigen::Matrix <value_type, 3, 1> vector3_t;
    } // namespace eigen
    typedef Eigen::Matrix <value_type, 5, 1> vector5_t;
    typedef Eigen::Matrix <value_type, 6, 1> vector6_t;
    typedef Eigen::Matrix <value_type, 7, 1> vector7_t;
    typedef Eigen::Quaternion<value_type> Quaternion_t;

    typedef pinocchio::ArrayXb ArrayXb;
    typedef ArrayXb bool_array_t;

    typedef std::pair<size_type, size_type> segment_t;
    typedef std::vector < segment_t > segments_t;

    HPP_PREDEF_CLASS (DistanceBetweenBodies);
    HPP_PREDEF_CLASS (DistanceBetweenPointsInBodies);
    HPP_PREDEF_CLASS (RelativeCom);
    HPP_PREDEF_CLASS (ComBetweenFeet);
    HPP_PREDEF_CLASS (StaticStability);
    HPP_PREDEF_CLASS (QPStaticStability);
    class ConvexShape;
    typedef std::vector <ConvexShape> ConvexShapes_t;
    HPP_PREDEF_CLASS (ConvexShapeContact);
    HPP_PREDEF_CLASS (ConvexShapeContactComplement);
    HPP_PREDEF_CLASS (ConvexShapeContactHold);
    HPP_PREDEF_CLASS (ConfigurationConstraint);
    HPP_PREDEF_CLASS (Identity);
    HPP_PREDEF_CLASS (AffineFunction);
    HPP_PREDEF_CLASS (ConstantFunction);

    typedef pinocchio::ObjectVector_t ObjectVector_t;
    typedef pinocchio::CollisionObjectPtr_t CollisionObjectPtr_t;
    typedef pinocchio::CollisionObjectConstPtr_t CollisionObjectConstPtr_t;
    typedef pinocchio::Configuration_t Configuration_t;
    typedef pinocchio::ConfigurationIn_t ConfigurationIn_t;
    typedef pinocchio::ConfigurationOut_t ConfigurationOut_t;
    typedef pinocchio::Device Device;
    typedef pinocchio::DevicePtr_t DevicePtr_t;
    typedef pinocchio::DeviceConstPtr_t DeviceConstPtr_t;
    typedef pinocchio::CenterOfMassComputation CenterOfMassComputation;
    typedef pinocchio::CenterOfMassComputationPtr_t CenterOfMassComputationPtr_t;
    typedef shared_ptr <DifferentiableFunction>
    DifferentiableFunctionPtr_t;
    typedef shared_ptr <DifferentiableFunctionSet>
    DifferentiableFunctionSetPtr_t;
    typedef DifferentiableFunctionSet DifferentiableFunctionStack
    HPP_CONSTRAINTS_DEPRECATED;
    typedef shared_ptr <ActiveSetDifferentiableFunction>
    ActiveSetDifferentiableFunctionPtr_t;
    typedef shared_ptr <DistanceBetweenBodies>
    DistanceBetweenBodiesPtr_t;
    typedef shared_ptr <DistanceBetweenPointsInBodies>
    DistanceBetweenPointsInBodiesPtr_t;
    typedef shared_ptr<RelativeCom> RelativeComPtr_t;
    typedef shared_ptr<ComBetweenFeet> ComBetweenFeetPtr_t;
    /// Plane polygon represented by its vertices
    /// Used to model contact surfaces for manipulation applications
    typedef std::vector<vector3_t> Shape_t;
    typedef std::pair <JointPtr_t, Shape_t> JointAndShape_t;
    typedef std::vector <JointAndShape_t> JointAndShapes_t;
    typedef shared_ptr<ConvexShapeContact>
      ConvexShapeContactPtr_t;
    typedef shared_ptr<ConvexShapeContactComplement>
      ConvexShapeContactComplementPtr_t;
    typedef shared_ptr<ConvexShapeContactHold>
      ConvexShapeContactHoldPtr_t;
    typedef shared_ptr<StaticStability> StaticStabilityPtr_t;
    typedef shared_ptr<QPStaticStability> QPStaticStabilityPtr_t;
    typedef shared_ptr<ConfigurationConstraint>
      ConfigurationConstraintPtr_t;
    typedef shared_ptr<Identity> IdentityPtr_t;
    typedef shared_ptr<AffineFunction> AffineFunctionPtr_t;
    typedef shared_ptr<ConstantFunction> ConstantFunctionPtr_t;

    template <int _Options> class GenericTransformation;

    /// \cond DEVEL
    const int RelativeBit       = 0x1;
    const int PositionBit       = 0x2;
    const int OrientationBit    = 0x4;
    const int OutputR3xSO3Bit      = 0x8;
    /// \endcond DEVEL
    typedef GenericTransformation<               PositionBit | OrientationBit > Transformation;
    typedef GenericTransformation<               PositionBit                  > Position;
    typedef GenericTransformation<                             OrientationBit > Orientation;
    typedef GenericTransformation< RelativeBit | PositionBit | OrientationBit > RelativeTransformation;
    typedef GenericTransformation< RelativeBit | PositionBit                  > RelativePosition;
    typedef GenericTransformation< RelativeBit |               OrientationBit > RelativeOrientation;
    typedef GenericTransformation<               PositionBit | OrientationBit | OutputR3xSO3Bit > TransformationR3xSO3;
    typedef GenericTransformation< RelativeBit | PositionBit | OrientationBit | OutputR3xSO3Bit > RelativeTransformationR3xSO3;
    typedef GenericTransformation<                             OrientationBit | OutputR3xSO3Bit > OrientationSO3;
    typedef GenericTransformation< RelativeBit |               OrientationBit | OutputR3xSO3Bit > RelativeOrientationSO3;

    typedef shared_ptr<Position> PositionPtr_t;
    typedef shared_ptr<Orientation> OrientationPtr_t;
    typedef shared_ptr<Transformation> TransformationPtr_t;
    typedef shared_ptr<RelativePosition> RelativePositionPtr_t;
    typedef shared_ptr<RelativeOrientation> RelativeOrientationPtr_t;
    typedef shared_ptr<RelativeTransformation>
      RelativeTransformationPtr_t;

    typedef Eigen::BlockIndex BlockIndex;

    HPP_PREDEF_CLASS (Implicit);
    typedef shared_ptr <Implicit> ImplicitPtr_t;
    typedef shared_ptr <const Implicit> ImplicitConstPtr_t;
    typedef std::vector < constraints::ImplicitPtr_t > NumericalConstraints_t;
    HPP_PREDEF_CLASS (ImplicitConstraintSet);
    typedef shared_ptr <ImplicitConstraintSet>
    ImplicitConstraintSetPtr_t;

    enum ComparisonType {
      Equality,
      EqualToZero,
      Superior,
      Inferior
    };
    typedef std::vector<ComparisonType> ComparisonTypes_t;

    HPP_PREDEF_CLASS (Explicit);
    typedef shared_ptr <Explicit> ExplicitPtr_t;
    typedef shared_ptr <const Explicit> ExplicitConstPtr_t;

    class ExplicitConstraintSet;
    namespace solver {
      class HierarchicalIterative;
      class BySubstitution;
    } // namespace solver

    namespace explicit_ {
      class Function;
      HPP_PREDEF_CLASS (RelativePose);
      typedef shared_ptr <RelativePose> RelativePosePtr_t;
      HPP_PREDEF_CLASS (RelativeTransformation);
      typedef shared_ptr <RelativeTransformation>
      RelativeTransformationPtr_t;
      HPP_PREDEF_CLASS (ImplicitFunction);
      typedef shared_ptr <ImplicitFunction> ImplicitFunctionPtr_t;
      HPP_PREDEF_CLASS(ConvexShapeContact);
      typedef shared_ptr<ConvexShapeContact> ConvexShapeContactPtr_t;
    } // namespace explicit_

    HPP_PREDEF_CLASS (LockedJoint);
    typedef shared_ptr <LockedJoint> LockedJointPtr_t;
    typedef shared_ptr <const LockedJoint> LockedJointConstPtr_t;
    typedef std::vector <LockedJointPtr_t> LockedJoints_t;

    namespace function {
      HPP_PREDEF_CLASS (OfParameterSubset);
      typedef shared_ptr <OfParameterSubset> OfParameterSubsetPtr_t;
    } // namespace function
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_FWD_HH
