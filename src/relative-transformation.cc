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

#include <hpp/model/device.hh>
#include <hpp/constraints/relative-transformation.hh>

namespace hpp {
  namespace constraints {
    RelativeTransformationPtr_t RelativeTransformation::create
    (const DevicePtr_t& robot, const JointPtr_t& joint1,
     const JointPtr_t& joint2, const Transform3f& reference)
    {
      RelativeTransformation* ptr = new RelativeTransformation
	(robot, joint1,	joint2, reference);
      RelativeTransformationPtr_t shPtr (ptr);
      return shPtr;
    }

    RelativeTransformation::RelativeTransformation
    (const DevicePtr_t& robot, const JointPtr_t& joint1,
     const JointPtr_t& joint2,
     const Transform3f& reference) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (), 6,
			      "RelativeTransformation"),
      relativeOrientation_ (robot, joint1, joint2, reference.getRotation (),
			    false),
      relativePosition_ (robot, joint1, joint2, reference.getTranslation (),
			 vector3_t (0, 0, 0)), reference_ (reference)
    {
    }

    void RelativeTransformation::impl_compute (vectorOut_t result,
					       ConfigurationIn_t argument)
      const throw ()
    {
      relativePosition_ (result.segment (0, 3), argument);
      relativeOrientation_ (result.segment (3, 3), argument);
    }

    void RelativeTransformation::impl_jacobian
    (matrixOut_t jacobian, ConfigurationIn_t arg) const throw ()
    {
      relativePosition_.jacobian (jacobian.topRows <3> (), arg);
      relativeOrientation_.jacobian (jacobian.bottomRows <3> (), arg);
    }

  } // namespace constraints
} // namespace hpp
