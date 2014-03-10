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

# include <hpp/core/fwd.hh>

namespace hpp {
  namespace constraints {
    HPP_PREDEF_CLASS (Orientation);
    HPP_PREDEF_CLASS (RelativeCom);
    HPP_PREDEF_CLASS (RelativeOrientation);
    HPP_PREDEF_CLASS (RelativePosition);
    HPP_PREDEF_CLASS (Position);

    typedef core::DevicePtr_t DevicePtr_t;
    typedef core::DifferentiableFunction DifferentiableFunction_t;
    typedef core::DifferentiableFunctionPtr_t DifferentiableFunctionPtr_t;
    typedef core::size_type size_type;
    typedef model::JointPtr_t JointPtr_t;
    typedef model::vector3_t vector3_t;
    typedef model::matrix3_t matrix3_t;
    typedef model::matrix_t matrix_t;
    typedef model::vector_t vector_t;
    typedef model::ComJacobian_t ComJacobian_t;
    typedef model::JointJacobian_t JointJacobian_t;
    typedef model::Transform3f Transform3f;
    typedef OrientationShPtr OrientationPtr_t;
    typedef PositionShPtr PositionPtr_t;
    typedef RelativeOrientationShPtr RelativeOrientationPtr_t;
    typedef RelativeComShPtr RelativeComPtr_t;
    typedef RelativePositionShPtr RelativePositionPtr_t;
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_FWD_HH
