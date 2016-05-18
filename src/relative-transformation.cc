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
#include <hpp/model/fcl-to-eigen.hh>
#include <hpp/model/joint.hh>
#include <hpp/constraints/tools.hh>
#include <hpp/constraints/relative-transformation.hh>
#include <hpp/constraints/macros.hh>

namespace hpp {
  namespace constraints {
    namespace deprecated {
    //using fcl::transpose;

    static size_type size (std::vector<bool> mask)
    {
      size_type res = 0;
      for (std::vector<bool>::iterator it = mask.begin (); it != mask.end ();
	   ++it) {
	if (*it) ++res;
      }
      return res;
    }

    static void cross (const vector3_t& v, eigen::matrix3_t& m)
    {
      m (0,1) = -v [2]; m (1,0) = v [2];
      m (0,2) = v [1]; m (2,0) = -v [1];
      m (1,2) = -v [0]; m (2,1) = v [0];
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

    RelativeTransformationPtr_t RelativeTransformation::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const JointPtr_t& joint2,
     const Transform3f& frame1, const Transform3f& frame2,
     std::vector <bool> mask)
    {
      RelativeTransformation* ptr = new RelativeTransformation
	(name, robot, joint1, joint2, frame1, frame2, mask);
      RelativeTransformationPtr_t shPtr (ptr);
      return shPtr;
    }

    RelativeTransformation::RelativeTransformation
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const JointPtr_t& joint2,
     const Transform3f& reference, std::vector <bool> mask) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
			      size (mask), name),
      robot_ (robot), joint1_ (joint1), joint2_ (joint2), F1inJ1_ (reference),
      F2inJ2_ (), mask_ (mask)
    {
      F2inJ2_.setIdentity ();
      value_.resize (6);
      jacobian_.resize (6, robot->numberDof ());
      cross1_.setZero ();
      cross2_.setZero ();
    }

    RelativeTransformation::RelativeTransformation
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const JointPtr_t& joint2,
     const Transform3f& frame1, const Transform3f& frame2,
     std::vector <bool> mask) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
			      size (mask), name),
      robot_ (robot), joint1_ (joint1), joint2_ (joint2), F1inJ1_ (frame1),
      F2inJ2_ (frame2), mask_ (mask)
    {
      value_.resize (6);
      jacobian_.resize (6, robot->numberDof ());
      cross1_.setZero ();
      cross2_.setZero ();
    }

    void RelativeTransformation::computeError (ConfigurationIn_t argument) const
    {
      hppDnum (info, "argument=" << argument.transpose ());
      if (argument.size () != latestArgument_.size () ||
	  argument != latestArgument_) {
	robot_->currentConfiguration (argument);
	robot_->computeForwardKinematics ();
	Transform3f J1;
	if (joint1_) {
	  J1 = joint1_->currentTransformation ();
	}
	const Transform3f& J2 = joint2_->currentTransformation ();
	Transform3f M1 (J1 * F1inJ1_);
	Transform3f M2 (J2 * F2inJ2_);
	hppDnum (info, "J1=" << J1);
	hppDnum (info, "J2=" << J2);
	hppDnum (info, "M1=" << M1);
	hppDnum (info, "M2=" << M2);
	Transform3f M (inverse (M1) * M2);
	const matrix3_t& Rerror (M.getRotation ());
	model::toEigen (M.getTranslation (), value_.head <3> ());
	hppDnum (info, "Rerror=" << Rerror);
	double tr = Rerror (0, 0) + Rerror (1, 1) + Rerror (2, 2);
	if (tr > 3) tr = 3;
	if (tr < -1) tr = -1;
	theta_ = acos ((tr - 1)/2);
	hppDnum (info, "theta_=" << theta_);
	assert (theta_ == theta_);
	if (theta_ < 1e-6) {
	  value_ [3] = (Rerror (2, 1) - Rerror (1, 2))/2;
	  value_ [4] = (Rerror (0, 2) - Rerror (2, 0))/2;
	  value_ [5] = (Rerror (1, 0) - Rerror (0, 1))/2;
	} else if (theta_ > M_PI - 1e-3) {
          const value_type phi = theta_ - M_PI;
          const value_type phi2_2 = phi * phi / 2;
          const value_type alpha = 1 - phi2_2;
          const value_type beta  = theta_*theta_ / ( 2 - phi2_2 );
          const value_type tmp0 = Rerror (0, 0) + alpha;
          const value_type tmp1 = Rerror (1, 1) + alpha;
          const value_type tmp2 = Rerror (2, 2) + alpha;
	  value_ [3] = (Rerror (2, 1) > Rerror (1, 2) ? 1 : -1 ) * (tmp0 > 0 ? sqrt(tmp0 * beta) : 0);
	  value_ [4] = (Rerror (0, 2) > Rerror (2, 0) ? 1 : -1 ) * (tmp1 > 0 ? sqrt(tmp1 * beta) : 0);
	  value_ [5] = (Rerror (1, 0) > Rerror (0, 1) ? 1 : -1 ) * (tmp2 > 0 ? sqrt(tmp2 * beta) : 0);
        } else {
	  value_ [3] = theta_*(Rerror (2, 1) -
			      Rerror (1, 2))/(2*sin(theta_));
	  value_ [4] = theta_*(Rerror (0, 2) -
			      Rerror (2, 0))/(2*sin(theta_));
	  value_ [5] = theta_*(Rerror (1, 0) -
			      Rerror (0, 1))/(2*sin(theta_));
	}
	latestArgument_ = argument;
      }
    }

    void RelativeTransformation::impl_compute (vectorOut_t result,
					       ConfigurationIn_t argument)
      const throw ()
    {
      computeError (argument);
      size_type index=0;
      for (size_type i=0; i<6; ++i) {
	if (mask_ [i]) {
	  result [index] = value_ [i]; ++index;
	}
      }
    }

    void RelativeTransformation::impl_jacobian
    (matrixOut_t jacobian, ConfigurationIn_t arg) const throw ()
    {
      computeError (arg);
      Transform3f J1;
      if (joint1_) {
	J1 = joint1_->currentTransformation ();
      }
      const Transform3f& J2 = joint2_->currentTransformation ();
      const vector3_t& t2inJ2 (F2inJ2_.getTranslation ());
      const matrix3_t& R1inJ1 (F1inJ1_.getRotation ());
      const vector3_t& t1 (J1.getTranslation ());
      const vector3_t& t2 (J2.getTranslation ());
      const matrix3_t& R1 (J1.getRotation ());
      const matrix3_t& R2 (J2.getRotation ());

      computeJlog (theta_, value_.tail <3> (), Jlog_);
      hppDnum (info, "Jlog_: " << Jlog_);
      //cross (R2*t2inJ2+t2-t1, cross1_);
      cross (R2*t2inJ2 + t2 - t1, cross1_);
      cross (R2*t2inJ2, cross2_);
      size_type leftCols = joint2_->jacobian().cols();
      jacobian_.rightCols (jacobian_.cols() - leftCols).setZero();
      if (joint1_) {
        jacobian_.topRows <3> ().leftCols (leftCols) =
#ifdef FCL_HAVE_EIGEN
          transpose (R1inJ1) * transpose (R1)
#else
          // This is a bug in Eigen that is fixed in a future version.
          // See https://bitbucket.org/eigen/eigen/commits/de7b8c9b1e86/
          transpose (R1*R1inJ1)
#endif
          *(cross1_*joint1_->jacobian ().bottomRows <3> ()-
				 cross2_*joint2_->jacobian ().bottomRows <3> ()+
				 joint2_->jacobian ().topRows <3>()-
				 joint1_->jacobian ().topRows <3>());
	jacobian_.bottomRows <3> ().leftCols (leftCols) =
	  Jlog_ *
#ifdef FCL_HAVE_EIGEN
          transpose (R1inJ1) * transpose (R1)
#else
          // This is a bug in Eigen that is fixed in a future version.
          // See https://bitbucket.org/eigen/eigen/commits/de7b8c9b1e86/
          transpose (R1*R1inJ1)
#endif
          *
	  (joint2_->jacobian ().bottomRows <3> () -
	   joint1_->jacobian ().bottomRows <3> ());
      } else {
	jacobian_.topRows <3> ().leftCols (leftCols) =
	  transpose (R1inJ1)*(-cross2_*joint2_->jacobian ().bottomRows <3> ()
			      + joint2_->jacobian ().topRows <3>());
	jacobian_.bottomRows <3> ().leftCols (leftCols) =
	  Jlog_ * transpose (R1inJ1) * (joint2_->jacobian ().bottomRows <3> ());
      }
      size_type index=0;
      for (size_type i=0; i<6; ++i) {
	if (mask_ [i]) {
	  jacobian.row (index) = jacobian_.row (i); ++index;
	}
      }
    }

    } // namespace deprecated
  } // namespace constraints
} // namespace hpp
