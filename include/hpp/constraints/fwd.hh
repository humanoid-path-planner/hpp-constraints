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

# include <hpp/model/fwd.hh>
# include <hpp/constraints/deprecated.hh>

namespace hpp {
  namespace constraints {
    HPP_PREDEF_CLASS (DifferentiableFunction);
    HPP_PREDEF_CLASS (DifferentiableFunctionStack);
    typedef model::size_type size_type;
    typedef model::value_type value_type;
    typedef model::JointPtr_t JointPtr_t;
    typedef model::vector3_t vector3_t;
    typedef model::matrix3_t matrix3_t;
    typedef model::matrix_t matrix_t;
    typedef Eigen::Ref <const matrix_t> matrixIn_t;
    typedef Eigen::Ref <matrix_t> matrixOut_t;
    typedef model::vector_t vector_t;
    typedef model::vectorIn_t vectorIn_t;
    typedef model::vectorOut_t vectorOut_t;
    typedef model::ComJacobian_t ComJacobian_t;
    typedef model::JointJacobian_t JointJacobian_t;
    typedef model::Transform3f Transform3f;
    namespace eigen {
      typedef Eigen::Matrix <value_type, 3, 3> matrix3_t;
      typedef Eigen::Matrix <value_type, 3, 1> vector3_t;
    } // namespace eigen
    typedef Eigen::Matrix <value_type, 5, 1> vector5_t;

    HPP_PREDEF_CLASS (DistanceBetweenBodies);
    HPP_PREDEF_CLASS (DistanceBetweenPointsInBodies);
    HPP_PREDEF_CLASS (RelativeCom);
    HPP_PREDEF_CLASS (ComBetweenFeet);
    HPP_PREDEF_CLASS (StaticStability);
    HPP_PREDEF_CLASS (QPStaticStability);
    HPP_PREDEF_CLASS (ConvexShapeContact);
    HPP_PREDEF_CLASS (ConvexShapeContactComplement);
    HPP_PREDEF_CLASS (ConfigurationConstraint);

    typedef model::ObjectVector_t ObjectVector_t;
    typedef model::CollisionObjectPtr_t CollisionObjectPtr_t;
    typedef model::Configuration_t Configuration_t;
    typedef model::ConfigurationIn_t ConfigurationIn_t;
    typedef model::ConfigurationOut_t ConfigurationOut_t;
    typedef model::Device Device;
    typedef model::DevicePtr_t DevicePtr_t;
    typedef model::CenterOfMassComputation CenterOfMassComputation;
    typedef model::CenterOfMassComputationPtr_t CenterOfMassComputationPtr_t;
    typedef boost::shared_ptr <DifferentiableFunction>
    DifferentiableFunctionPtr_t;
    typedef boost::shared_ptr <DifferentiableFunctionStack>
    DifferentiableFunctionStackPtr_t;
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
    /// \endcond DEVEL
    typedef GenericTransformation<               PositionBit | OrientationBit > Transformation;
    typedef GenericTransformation<               PositionBit                  > Position;
    typedef GenericTransformation<                             OrientationBit > Orientation;
    typedef GenericTransformation< RelativeBit | PositionBit | OrientationBit > RelativeTransformation;
    typedef GenericTransformation< RelativeBit | PositionBit                  > RelativePosition;
    typedef GenericTransformation< RelativeBit |               OrientationBit > RelativeOrientation;
    typedef boost::shared_ptr<Position> PositionPtr_t;
    typedef boost::shared_ptr<Orientation> OrientationPtr_t;
    typedef boost::shared_ptr<Transformation> TransformationPtr_t;
    typedef boost::shared_ptr<RelativePosition> RelativePositionPtr_t;
    typedef boost::shared_ptr<RelativeOrientation> RelativeOrientationPtr_t;
    typedef boost::shared_ptr<RelativeTransformation>
      RelativeTransformationPtr_t;

    namespace deprecated {
      HPP_PREDEF_CLASS (Position);
      HPP_PREDEF_CLASS (Orientation);
      HPP_PREDEF_CLASS (Transformation);
      HPP_PREDEF_CLASS (RelativeOrientation);
      HPP_PREDEF_CLASS (RelativePosition);
      HPP_PREDEF_CLASS (RelativeTransformation);
      typedef boost::shared_ptr<Position> PositionPtr_t;
      typedef boost::shared_ptr<Orientation> OrientationPtr_t;
      typedef boost::shared_ptr<Transformation> TransformationPtr_t;
      typedef boost::shared_ptr<RelativePosition> RelativePositionPtr_t;
      typedef boost::shared_ptr<RelativeOrientation> RelativeOrientationPtr_t;
      typedef boost::shared_ptr<RelativeTransformation>
        RelativeTransformationPtr_t;
    } // namespace deprecated
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_FWD_HH
