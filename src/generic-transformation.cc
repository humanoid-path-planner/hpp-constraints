// Copyright (c) 2016, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//
// This file is part of hpp-constraints.
// hpp-constraints is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-constraints is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-constraints. If not, see <http://www.gnu.org/licenses/>.

#include <hpp/constraints/generic-transformation.hh>

#include <boost/math/constants/constants.hpp>

#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>

#include <hpp/constraints/macros.hh>

namespace hpp {
  namespace constraints {
    namespace {
      static inline size_type size (std::vector<bool> mask)
      {
        size_type res = 0;
        for (std::vector<bool>::iterator it = mask.begin (); it != mask.end ();
            ++it) {
          if (*it) ++res;
        }
        return res;
      }

      template <bool flag> struct unary
      {
        template <bool rel, bool pos> static inline void log (
            const Transform3f&, const GenericTransformationData<rel, pos, flag>&) {}
        template <bool rel, bool pos> static inline void Jlog (
            const GenericTransformationData<rel, pos, flag>&) {}
      };
      template <> struct unary <true>
      {
        template <bool rel, bool pos> static inline void log (
            const Transform3f& M, const GenericTransformationData<rel, pos, true>& d)
          {
            const matrix3_t& Rerror (M.getRotation ());
            value_type tr = Rerror.trace();
            if (tr > 3)       d.theta = 0; // acos((3-1)/2)
            else if (tr < -1) d.theta = ::boost::math::constants::pi<value_type>(); // acos((-1-1)/2)
            else              d.theta = acos ((tr - 1)/2);
            hppDnum (info, "theta_=" << d.theta);
            assert (d.theta == d.theta);
            const value_type t = ((d.theta > 1e-6)? d.theta / sin(d.theta) : 1) / 2;
            d.value((pos?3:0)+0) = t * (Rerror (2, 1) - Rerror (1, 2));
            d.value((pos?3:0)+1) = t * (Rerror (0, 2) - Rerror (2, 0));
            d.value((pos?3:0)+2) = t * (Rerror (1, 0) - Rerror (0, 1));
          }
        template <bool rel, bool pos> static inline void Jlog (
            const GenericTransformationData<rel, pos, true>& d)
          {
            if (d.theta < 1e-6) d.Jlog.setIdentity ();
            else {
              // Jlog = alpha I
              const value_type ct = cos(d.theta), st = sin(d.theta);
              const value_type st_1mct = st/(1-ct);

              d.Jlog.setZero ();
              d.Jlog.diagonal().setConstant (d.theta*st_1mct);

              // Jlog += -r_{\times}/2
              d.Jlog(0,1) =  d.value((pos?3:0)+2); d.Jlog(1,0) = -d.value((pos?3:0)+2);
              d.Jlog(0,2) = -d.value((pos?3:0)+1); d.Jlog(2,0) =  d.value((pos?3:0)+1);
              d.Jlog(1,2) =  d.value((pos?3:0)+0); d.Jlog(2,1) = -d.value((pos?3:0)+0);
              d.Jlog /= 2;

              const value_type alpha = 1/(d.theta*d.theta) - st_1mct/(2*d.theta);
              d.Jlog.noalias() += alpha * d.value.template tail<3>() * d.value.template tail<3>().transpose ();
            }
          } 
      };

      template <bool lflag /*rel*/, bool rflag /*false*/> struct binary
      {
        // the first template allow us to consider relative transformation as
        // absolute when joint1 is NULL, at run time
        template <bool rel, bool pos> static inline void Jorientation (
            const GenericTransformationData<rel, pos, rflag>&) {}
        template <bool rel, bool ori> static inline void Jtranslation (
            const GenericTransformationData<rel, rflag, ori>&) {}
      };
      template <> struct binary<false, true> // Absolute
      {
        template <bool rel, bool pos> static inline void Jorientation (
            const GenericTransformationData<rel, pos, true>& d)
        {
          const matrix3_t& R1inJ1 (d.F1inJ1.getRotation ());
          const size_type leftCols = d.joint2->jacobian().cols();
          d.jacobian.template bottomRows<3>().leftCols (leftCols) =
            d.Jlog * transpose(R1inJ1) * (d.joint2->jacobian().template bottomRows<3>());
        }
        template <bool rel, bool ori> static inline void Jtranslation (
            const GenericTransformationData<rel, true, ori>& d)
        {
          const matrix3_t& R1inJ1 (d.F1inJ1.getRotation ());
          const size_type leftCols = d.joint2->jacobian().cols();
          d.jacobian.template topRows<3>().leftCols (leftCols) =
              transpose(R1inJ1) * (
                    d.joint2->jacobian().template bottomRows<3>().colwise().cross(d.cross2)
                  + d.joint2->jacobian().template topRows<3>());
        }
      };
      template <> struct binary<true, true> // Relative
      {
        template <bool pos> static inline void Jorientation (
            const GenericTransformationData<true, pos, true>& d)
        {
          const matrix3_t& R1inJ1 (d.F1inJ1.getRotation ());
          const Transform3f& J1 = d.joint1->currentTransformation ();
          const matrix3_t& R1 (J1.getRotation ());
          const size_type leftCols = d.joint2->jacobian().cols();
          d.jacobian.template bottomRows<3>().leftCols (leftCols) =
            d.Jlog * transpose (R1inJ1) * transpose (R1) *
            (  d.joint2->jacobian().template bottomRows<3>()
             - d.joint1->jacobian().template bottomRows<3>());
        }
        template <bool ori> static inline void Jtranslation (
            const GenericTransformationData<true, true, ori>& d)
        {
          const matrix3_t& R1inJ1 (d.F1inJ1.getRotation ());
          const Transform3f& J1 = d.joint1->currentTransformation ();
          const matrix3_t& R1 (J1.getRotation ());
          const size_type leftCols = d.joint2->jacobian().cols();
          // This is a bug in Eigen that is fixed in a future version.
          // See https://bitbucket.org/eigen/eigen/commits/de7b8c9b1e86/
          d.jacobian.template topRows<3>().leftCols (leftCols) =
            transpose(R1inJ1) * transpose(R1) * (
                - d.joint1->jacobian().template bottomRows<3>().colwise().cross(d.cross1)
                + d.joint2->jacobian().template bottomRows<3>().colwise().cross(d.cross2)
                + d.joint2->jacobian().template topRows<3>()
                - d.joint1->jacobian().template topRows<3>());
        }
      };

      template <bool rel> static inline void relativeTransform (
            const JointPtr_t j1, const Transform3f& f1,
            const JointPtr_t j2, const Transform3f& f2,
            Transform3f& out)
        {
          // TODO: when position only, do not compute the relative orientation.
          const Transform3f& J2 = j2->currentTransformation ();
          if (rel && j1) {
            const Transform3f& J1 = j1->currentTransformation ();
            out = f1.inverseTimes(J1.inverseTimes(J2 * f2));
          } else {
            out = f1.inverseTimes(J2 * f2);
          }
        }

      template <bool rel, bool pos, bool ori> struct compute
      {
        static inline void error (const GenericTransformationData<rel, pos, ori>& d)
        {
          Transform3f M;
          relativeTransform<rel> (d.getJoint1(), d.F1inJ1, d.joint2, d.F2inJ2, M);
          if (pos) d.value.template head<3>() = M.getTranslation();
          unary<ori>::log(M, d);
        }

        static inline void jacobian (const GenericTransformationData<rel, pos, ori>& d)
        {
          const Transform3f& J2 = d.joint2->currentTransformation ();
          const vector3_t& t2inJ2 (d.F2inJ2.getTranslation ());
          const vector3_t& t2 (J2.getTranslation ());
          const matrix3_t& R2 (J2.getRotation ());

          d.cross2.noalias() = R2*t2inJ2;

          unary<ori>::Jlog (d);
          hppDnum (info, "Jlog_: " << d.Jlog);

          const size_type leftCols = d.joint2->jacobian().cols();
          d.jacobian.rightCols (d.jacobian.cols() - leftCols).setZero();

          // rel:           relative known at compile time
          // d.getJoint1(): relative known at run time
          if (rel && d.getJoint1()) {
            const Transform3f& J1 = d.getJoint1()->currentTransformation ();
            const vector3_t& t1 (J1.getTranslation ());
            d.cross1.noalias() = d.cross2 + t2 - t1;
            binary<rel, pos>::Jtranslation (d);
            binary<rel, ori>::Jorientation (d);
          } else {
            d.cross1.noalias() = d.cross2 + t2;
            binary<false, pos>::Jtranslation (d);
            binary<false, ori>::Jorientation (d);
          }
        }
      };
    }

    template <int _Options> typename GenericTransformation<_Options>::Ptr_t
      GenericTransformation<_Options>::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint2,
     const Transform3f& reference, std::vector <bool> mask)
    {
      GenericTransformation<_Options>* ptr =
        new GenericTransformation<_Options> (name, robot, mask);
      ptr->joint1 (NULL);
      ptr->joint2 (joint2);
      ptr->reference (reference);
      Ptr_t shPtr (ptr);
      ptr->init (shPtr);
      return shPtr;
    }

    template <int _Options> typename GenericTransformation<_Options>::Ptr_t
      GenericTransformation<_Options>::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint2,
     const Transform3f& frame1, const Transform3f& frame2,
     std::vector <bool> mask)
    {
      GenericTransformation<_Options>* ptr =
        new GenericTransformation<_Options> (name, robot, mask);
      ptr->joint1 (NULL);
      ptr->joint2 (joint2);
      ptr->frame1InJoint1 (frame1);
      ptr->frame2InJoint2 (frame2);
      Ptr_t shPtr (ptr);
      return shPtr;
    }

    template <int _Options> typename GenericTransformation<_Options>::Ptr_t
      GenericTransformation<_Options>::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const JointPtr_t& joint2,
     const Transform3f& reference, std::vector <bool> mask)
    {
      GenericTransformation<_Options>* ptr =
        new GenericTransformation<_Options> (name, robot, mask);
      ptr->joint1 (joint1);
      ptr->joint2 (joint2);
      ptr->reference (reference);
      Ptr_t shPtr (ptr);
      ptr->init (shPtr);
      return shPtr;
    }

    template <int _Options> typename GenericTransformation<_Options>::Ptr_t
      GenericTransformation<_Options>::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointPtr_t& joint1, const JointPtr_t& joint2,
     const Transform3f& frame1, const Transform3f& frame2,
     std::vector <bool> mask)
    {
      GenericTransformation<_Options>* ptr =
        new GenericTransformation<_Options> (name, robot, mask);
      ptr->joint1 (joint1);
      ptr->joint2 (joint2);
      ptr->frame1InJoint1 (frame1);
      ptr->frame2InJoint2 (frame2);
      Ptr_t shPtr (ptr);
      return shPtr;
    }

    template <int _Options>
      GenericTransformation<_Options>::GenericTransformation
      (const std::string& name, const DevicePtr_t& robot,
       std::vector <bool> mask) :
        DifferentiableFunction (robot->configSize (), robot->numberDof (),
			      size (mask), name),
      robot_ (robot), d_(robot->numberDof()), mask_ (mask)
    {
      assert(mask.size()==ValueSize);
    }

    template <int _Options>
    void GenericTransformation<_Options>::computeError (ConfigurationIn_t argument) const
    {
      hppDnum (info, "argument=" << argument.transpose ());
      if (argument.size () != latestArgument_.size () ||
	  argument != latestArgument_) {
	robot_->currentConfiguration (argument);
	robot_->computeForwardKinematics ();
        compute<IsRelative, ComputePosition, ComputeOrientation>::error (d_);
	latestArgument_ = argument;
      }
    }

    template <int _Options>
    void GenericTransformation<_Options>::impl_compute (vectorOut_t result,
					       ConfigurationIn_t argument)
      const throw ()
    {
      computeError (argument);
      size_type index=0;
      for (size_type i=0; i<ValueSize; ++i) {
	if (mask_ [i]) {
	  result [index] = d_.value[i]; ++index;
	}
      }
    }

    template <int _Options>
    void GenericTransformation<_Options>::impl_jacobian
    (matrixOut_t jacobian, ConfigurationIn_t arg) const throw ()
    {
      computeError (arg);
      compute<IsRelative, ComputePosition, ComputeOrientation>::jacobian (d_);
      size_type index=0;
      for (size_type i=0; i<ValueSize; ++i) {
	if (mask_ [i]) {
	  jacobian.row (index) = d_.jacobian.row(i); ++index;
	}
      }
    }

    /// Force instanciation of relevant classes
    template class GenericTransformation<               PositionBit | OrientationBit >;
    template class GenericTransformation<               PositionBit                  >;
    template class GenericTransformation<                             OrientationBit >;
    template class GenericTransformation< RelativeBit | PositionBit | OrientationBit >;
    template class GenericTransformation< RelativeBit | PositionBit                  >;
    template class GenericTransformation< RelativeBit |               OrientationBit >;
  } // namespace constraints
} // namespace hpp
