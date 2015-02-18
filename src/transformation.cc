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

#include <Eigen/Dense>
#include <hpp/model/device.hh>
#include <hpp/constraints/orientation.hh>
#include <hpp/constraints/transformation.hh>

namespace hpp {
  namespace constraints {
    static matrix3_t identity_;
    static size_type size (std::vector<bool> mask)
    {
      size_type res = 0;
      for (std::vector<bool>::iterator it = mask.begin (); it != mask.end ();
	   ++it) {
	if (*it) ++res;
      }
      return res;
    }

    TransformationPtr_t Transformation::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1,
     const Transform3f& reference, std::vector <bool> mask)
    {
      identity_.setIdentity ();
      Transformation* ptr = new Transformation
	(name, robot, joint1, reference, mask);
      TransformationPtr_t shPtr (ptr);
      return shPtr;
    }

    TransformationPtr_t Transformation::create
    (const DevicePtr_t& robot, const JointPtr_t& joint1,
     const Transform3f& reference, std::vector <bool> mask)
    {
      identity_.setIdentity ();
      Transformation* ptr = new Transformation
	("Transformation", robot, joint1, reference, mask);
      TransformationPtr_t shPtr (ptr);
      return shPtr;
    }

    Transformation::Transformation
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1,
     const Transform3f& reference, std::vector <bool> mask) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
			      size (mask), name),
      orientation_ ("", robot, joint1, reference.getRotation (),
		    boost::assign::list_of (mask [3])(mask [4])
		    (mask [5])),
      position_ ("", robot, joint1, vector3_t (0, 0, 0),
		 reference.getTranslation (), identity_,
		 boost::assign::list_of (mask [0])(mask [1])
		 (mask [2])), reference_ (reference)
    {
      sizeTranslation_ = position_.outputSize ();
      sizeOrientation_ = orientation_.outputSize ();
    }

    void Transformation::impl_compute (vectorOut_t result,
				       ConfigurationIn_t argument)
      const throw ()
    {
      position_ (result.segment (0, sizeTranslation_), argument);
      orientation_ (result.segment (sizeTranslation_,
					    sizeOrientation_), argument);
    }

    void Transformation::impl_jacobian
    (matrixOut_t jacobian, ConfigurationIn_t arg) const throw ()
    {
      position_.jacobian (jacobian.topRows (sizeTranslation_), arg);
      orientation_.jacobian (jacobian.bottomRows (sizeOrientation_), arg);
    }

  } // namespace constraints
} // namespace hpp
