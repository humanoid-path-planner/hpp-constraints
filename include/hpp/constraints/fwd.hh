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

#ifndef HPP_CONSTRAINTS_FWD_HH
# define HPP_CONSTRAINTS_FWD_HH

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
    typedef pinocchio::vector3_t vector3_t;
    typedef pinocchio::matrix3_t matrix3_t;
    typedef Eigen::Matrix<value_type, 6, 6> matrix6_t;
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
    HPP_PREDEF_CLASS (ConvexShapeContact);
    HPP_PREDEF_CLASS (ConvexShapeContactComplement);
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
    typedef pinocchio::CenterOfMassComputation CenterOfMassComputation;
    typedef pinocchio::CenterOfMassComputationPtr_t CenterOfMassComputationPtr_t;
    typedef boost::shared_ptr <DifferentiableFunction>
    DifferentiableFunctionPtr_t;
    typedef boost::shared_ptr <DifferentiableFunctionSet>
    DifferentiableFunctionSetPtr_t;
    typedef DifferentiableFunctionSet DifferentiableFunctionStack
    HPP_CONSTRAINTS_DEPRECATED;
    typedef boost::shared_ptr <ActiveSetDifferentiableFunction>
    ActiveSetDifferentiableFunctionPtr_t;
    typedef boost::shared_ptr <DistanceBetweenBodies>
    DistanceBetweenBodiesPtr_t;
    typedef boost::shared_ptr <DistanceBetweenPointsInBodies>
    DistanceBetweenPointsInBodiesPtr_t;
    typedef boost::shared_ptr<RelativeCom> RelativeComPtr_t;
    typedef boost::shared_ptr<ComBetweenFeet> ComBetweenFeetPtr_t;
    typedef boost::shared_ptr<ConvexShapeContact>
      ConvexShapeContactPtr_t;
    typedef boost::shared_ptr<ConvexShapeContactComplement>
      ConvexShapeContactComplementPtr_t;
    typedef boost::shared_ptr<StaticStability> StaticStabilityPtr_t;
    typedef boost::shared_ptr<QPStaticStability> QPStaticStabilityPtr_t;
    typedef boost::shared_ptr<ConfigurationConstraint>
      ConfigurationConstraintPtr_t;
    typedef boost::shared_ptr<Identity> IdentityPtr_t;
    typedef boost::shared_ptr<AffineFunction> AffineFunctionPtr_t;
    typedef boost::shared_ptr<ConstantFunction> ConstantFunctionPtr_t;

    typedef HPP_CONSTRAINTS_DEPRECATED ConvexShapeContact StaticStabilityGravity;
    typedef HPP_CONSTRAINTS_DEPRECATED ConvexShapeContactComplement StaticStabilityGravityComplement;
    typedef HPP_CONSTRAINTS_DEPRECATED
      ConvexShapeContactPtr_t StaticStabilityGravityPtr_t;
    typedef HPP_CONSTRAINTS_DEPRECATED
      ConvexShapeContactComplementPtr_t StaticStabilityGravityComplementPtr_t;

    template <int _Options> class GenericTransformation;

    /// \cond DEVEL
    const int RelativeBit       = 0x1;
    const int PositionBit       = 0x2;
    const int OrientationBit    = 0x4;
    const int OutputSE3Bit      = 0x8;
    /// \endcond DEVEL
    typedef GenericTransformation<               PositionBit | OrientationBit > Transformation;
    typedef GenericTransformation<               PositionBit                  > Position;
    typedef GenericTransformation<                             OrientationBit > Orientation;
    typedef GenericTransformation< RelativeBit | PositionBit | OrientationBit > RelativeTransformation;
    typedef GenericTransformation< RelativeBit | PositionBit                  > RelativePosition;
    typedef GenericTransformation< RelativeBit |               OrientationBit > RelativeOrientation;
    typedef GenericTransformation<               PositionBit | OrientationBit | OutputSE3Bit > TransformationSE3;
    typedef GenericTransformation< RelativeBit | PositionBit | OrientationBit | OutputSE3Bit > RelativeTransformationSE3;
    typedef GenericTransformation<                             OrientationBit | OutputSE3Bit > OrientationSO3;
    typedef GenericTransformation< RelativeBit |               OrientationBit | OutputSE3Bit > RelativeOrientationSO3;

    typedef boost::shared_ptr<Position> PositionPtr_t;
    typedef boost::shared_ptr<Orientation> OrientationPtr_t;
    typedef boost::shared_ptr<Transformation> TransformationPtr_t;
    typedef boost::shared_ptr<RelativePosition> RelativePositionPtr_t;
    typedef boost::shared_ptr<RelativeOrientation> RelativeOrientationPtr_t;
    typedef boost::shared_ptr<RelativeTransformation>
      RelativeTransformationPtr_t;

    typedef Eigen::BlockIndex BlockIndex;

    HPP_PREDEF_CLASS (Implicit);
    typedef boost::shared_ptr <Implicit> ImplicitPtr_t;
    typedef boost::shared_ptr <const Implicit> ImplicitConstPtr_t;
    typedef std::vector < constraints::ImplicitPtr_t > NumericalConstraints_t;
    HPP_PREDEF_CLASS (ImplicitConstraintSet);
    typedef boost::shared_ptr <ImplicitConstraintSet>
    ImplicitConstraintSetPtr_t;

    enum ComparisonType {
      Equality,
      EqualToZero,
      Superior,
      Inferior
    };
    typedef std::vector<ComparisonType> ComparisonTypes_t;

    HPP_PREDEF_CLASS (Explicit);
    typedef boost::shared_ptr <Explicit> ExplicitPtr_t;
    typedef boost::shared_ptr <const Explicit> ExplicitConstPtr_t;

    class ExplicitConstraintSet;
    namespace solver {
      class HierarchicalIterative;
      class BySubstitution;
    } // namespace solver

    namespace implicit {
      HPP_PREDEF_CLASS (RelativePose);
      typedef boost::shared_ptr <RelativePose> RelativePosePtr_t;
      typedef boost::shared_ptr <const RelativePose> RelativePoseConstPtr_t;
    } // namespace implicit

    namespace explicit_ {
      class Function;
      HPP_PREDEF_CLASS (RelativePose);
      typedef boost::shared_ptr <RelativePose> RelativePosePtr_t;
      HPP_PREDEF_CLASS (RelativeTransformation);
      typedef boost::shared_ptr <RelativeTransformation>
      RelativeTransformationPtr_t;
      HPP_PREDEF_CLASS (ImplicitFunction);
      typedef boost::shared_ptr <ImplicitFunction> ImplicitFunctionPtr_t;
    } // namespace explicit

    HPP_PREDEF_CLASS (LockedJoint);
    typedef boost::shared_ptr <LockedJoint> LockedJointPtr_t;
    typedef boost::shared_ptr <const LockedJoint> LockedJointConstPtr_t;
    typedef std::vector <LockedJointPtr_t> LockedJoints_t;

    namespace function {
      HPP_PREDEF_CLASS (OfParameterSubset);
      typedef boost::shared_ptr <OfParameterSubset> OfParameterSubsetPtr_t;
    } // namespace function
    typedef ExplicitConstraintSet ExplicitSolver HPP_CONSTRAINTS_DEPRECATED;
    typedef solver::HierarchicalIterative HierarchicalIterativeSolver
    HPP_CONSTRAINTS_DEPRECATED;
    typedef solver::BySubstitution HybridSolver HPP_CONSTRAINTS_DEPRECATED;
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_FWD_HH
