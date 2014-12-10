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
#include <hpp/constraints/orientation.hh>
#include <hpp/constraints/relative-transformation.hh>

namespace hpp {
  namespace constraints {
    static size_type size (std::vector<bool> mask)
    {
      size_type res = 0;
      for (std::vector<bool>::iterator it = mask.begin (); it != mask.end ();
	   ++it) {
	if (*it) ++res;
      }
      return res;
    }

    RelativeTransformationPtr_t RelativeTransformation::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const JointPtr_t& joint2,
     const Transform3f& reference, std::vector <bool> mask)
    {
      RelativeTransformation* ptr = new RelativeTransformation
	(name, robot, joint1, joint2, reference, mask);
      RelativeTransformationPtr_t shPtr (ptr);
      return shPtr;
    }

    RelativeTransformationPtr_t RelativeTransformation::create
    (const DevicePtr_t& robot, const JointPtr_t& joint1,
     const JointPtr_t& joint2, const Transform3f& reference,
     std::vector <bool> mask)
    {
      RelativeTransformation* ptr = new RelativeTransformation
	("RelativeTransformation", robot, joint1, joint2, reference, mask);
      RelativeTransformationPtr_t shPtr (ptr);
      return shPtr;
    }

    RelativeTransformation::RelativeTransformation
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const JointPtr_t& joint2,
     const Transform3f& reference, std::vector <bool> mask) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
			      size (mask), name),
      relativeOrientation_ ("", robot, joint1, joint2, reference.getRotation (),
			    boost::assign::list_of (mask [3])(mask [4])
			    (mask [5])),
      relativePosition_ ("", robot, joint1, joint2, reference.getTranslation (),
			 vector3_t (0, 0, 0),
			 boost::assign::list_of (mask [0])(mask [1])
			 (mask [2])), reference_ (reference)
    {
      sizeTranslation_ = relativePosition_.outputSize ();
      sizeOrientation_ = relativeOrientation_.outputSize ();
    }

    void RelativeTransformation::impl_compute (vectorOut_t result,
					       ConfigurationIn_t argument)
      const throw ()
    {
      relativePosition_ (result.segment (0, sizeTranslation_), argument);
      relativeOrientation_ (result.segment (sizeTranslation_,
					    sizeOrientation_), argument);
    }

    void RelativeTransformation::impl_jacobian
    (matrixOut_t jacobian, ConfigurationIn_t arg) const throw ()
    {
      relativePosition_.jacobian (jacobian.topRows (sizeTranslation_), arg);
      relativeOrientation_.jacobian (jacobian.bottomRows (sizeOrientation_),
				     arg);
    }

  } // namespace constraints
} // namespace hpp
