// Copyright (c) 2014, LAAS-CNRS
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

#ifndef HPP_CONSTRAINTS_TOOL_HH
#define HPP_CONSTRAINTS_TOOL_HH

#include <boost/math/constants/constants.hpp>

#include <pinocchio/spatial/se3.hpp>

#include <hpp/constraints/fwd.hh>

namespace hpp {
  namespace constraints {

    template < typename VectorType, typename MatrixType >
    static void computeCrossMatrix (const VectorType& v, MatrixType& m)
    {
      m.diagonal ().setZero ();
      m (0,1) = -v [2]; m (1,0) = v [2];
      m (0,2) = v [1]; m (2,0) = -v [1];
      m (1,2) = -v [0]; m (2,1) = v [0];
    }

    /// \addtogroup hpp_constraints_tools
    /// \{

    /// Compute log of rotation matrix as a 3d vector
    ///
    /// \param R rotation matrix in SO(3),
    /// \retval theta angle of rotation,
    /// \retval result 3d vector \f$\mathbf{r}\f$ such that
    ///         \f$R=\exp [\mathbf{r}]_{\times}\f$
    template <typename Derived> inline void logSO3
    (const matrix3_t& R, value_type& theta,
     Eigen::MatrixBase<Derived> const& result)
    {
      Eigen::MatrixBase<Derived>& value = const_cast<Eigen::MatrixBase<Derived>&> (result);
      const value_type PI = ::boost::math::constants::pi<value_type>();
      const value_type tr = R.trace();
      if (tr > 3)       theta = 0; // acos((3-1)/2)
      else if (tr < -1) theta = PI; // acos((-1-1)/2)
      else              theta = acos ((tr - 1)/2);
      assert (theta == theta);
      // From runs of tests/logarithm.cc: 1e-6 is too small.
      if (theta < PI - 1e-2) {
        const value_type t = ((theta > 1e-6)? theta / sin(theta) : 1) / 2;
        value(0) = t * (R (2, 1) - R (1, 2));
        value(1) = t * (R (0, 2) - R (2, 0));
        value(2) = t * (R (1, 0) - R (0, 1));
      } else {
        // 1e-2: A low value is not required since the computation is
        // using explicit formula. However, the precision of this method
        // is the square root of the precision with the antisymmetric
        // method (Nominal case).
        const value_type cphi = cos(theta - PI);
        const value_type beta  = theta*theta / ( 1 + cphi );
        const value_type tmp0 = (R (0, 0) + cphi) * beta;
        const value_type tmp1 = (R (1, 1) + cphi) * beta;
        const value_type tmp2 = (R (2, 2) + cphi) * beta;
        value(0) = (R (2, 1) > R (1, 2) ? 1 : -1) * (tmp0 > 0 ? sqrt(tmp0) : 0);
        value(1) = (R (0, 2) > R (2, 0) ? 1 : -1) * (tmp1 > 0 ? sqrt(tmp1) : 0);
        value(2) = (R (1, 0) > R (0, 1) ? 1 : -1) * (tmp2 > 0 ? sqrt(tmp2) : 0);
      }
    }


    /** Compute jacobian of function log of rotation matrix in SO(3)

        Let us consider a matrix
        \f$R=\exp \left[\mathbf{r}\right]_{\times}\in SO(3)\f$.
        This functions computes the Jacobian of the function from
        \f$SO(3)\f$ into \f$\mathbf{R}^3\f$ that maps \f$R\f$ to
        \f$\mathbf{r}\f$. In other words,
        \f{equation*}
        \dot{\mathbf{r}} = J_{log}(R)\ \omega\,\,\,\mbox{with}\,\,\,
        \dot {R} = \left[\omega\right]_{\times} R
        \f}
        \warning Two representations of the angular velocity \f$\omega\f$ are
                 possible:
                 \li \f$\dot{R} = \left[\omega\right]_{\times}R\f$ or
                 \li \f$\dot{R} = R\left[\omega\right]_{\times}\f$.

                 The expression below assumes the second representation is used.
        \param theta angle of rotation \f$R\f$, also \f$\|r\|\f$,
        \param log 3d vector \f$\mathbf{r}\f$,
        \retval Jlog matrix \f$J_{log} (R)\f$.

        \f{align*}
        J_{log} (R) &=& \frac{\|\mathbf{r}\|\sin\|\mathbf{r}\|}{2(1-\cos\|\mathbf{r}\|)} I_3 - \frac {1}{2}\left[\mathbf{r}\right]_{\times} + (\frac{1}{\|\mathbf{r}\|^2} - \frac{\sin\|\mathbf{r}\|}{2\|\mathbf{r}\|(1-\cos\|\mathbf{r}\|)}) \mathbf{r}\mathbf{r}^T\\
         &=& I_3 +\frac{1}{2}\left[\mathbf{r}\right]_{\times} +  \left(\frac{2(1-\cos\|\mathbf{r}\|) - \|\mathbf{r}\|\sin\|\mathbf{r}\|}{2\|\mathbf{r}\|^2(1-\cos\|\mathbf{r}\|)}\right)\left[\mathbf{r}\right]_{\times}^2
        \f}
    */
    template <typename Derived>
      void JlogSO3 (const value_type& theta,
                    const Eigen::MatrixBase<Derived>& log, matrix3_t& Jlog)
    {
      if (theta < 1e-6)
        Jlog.setIdentity();
      else {
        // Jlog = alpha I
        const value_type ct = cos(theta), st = sin(theta);
        const value_type st_1mct = st/(1-ct);

        Jlog.setZero ();
        Jlog.diagonal().setConstant (theta*st_1mct);

        // Jlog += r_{\times}/2
        Jlog(0,1) = -log(2); Jlog(1,0) =  log(2);
        Jlog(0,2) =  log(1); Jlog(2,0) = -log(1);
        Jlog(1,2) = -log(0); Jlog(2,1) =  log(0);
        Jlog /= 2;

        const value_type alpha = 1/(theta*theta) - st_1mct/(2*theta);
        Jlog.noalias() += alpha * log * log.transpose ();
      }
    }

    /// Compute log of rigid-body transform
    ///
    /// \param M rigid body transform,
    /// \retval theta angle of rotation,
    /// \retval result 6d vector \f$(\mathbf{v},\mathbf{r})\f$ such that
    /// the screw motion of linear velocity (of the origin) \f$\mathbf{v}\f$
    /// expressed in the moving frame and of angular velocity \f$\mathbf{r}\f$
    /// expressed in the moving reaches $M$ in unit time.
    template <typename Derived> inline void logSE3
    (const Transform3f& M, Eigen::MatrixBase<Derived>& result)
    {
      assert (result.size () == 6);
      Eigen::MatrixBase<Derived>& value =
        const_cast<Eigen::MatrixBase<Derived>&> (result);
      const matrix3_t & R = M.rotation();
      const vector3_t & p = M.translation();
      value_type theta;
      vector3_t r;
      logSO3 (R, theta, r);
      value.segment (3, 3) = r;
      value_type alpha, beta;
      if (fabs (theta) < 1e-2) {
        alpha = 1 - theta*theta/12 - theta*theta*theta*theta/720;
        beta = 1./12 + theta*theta/720;
      } else {
        alpha = theta*sin(theta)/(2*(1-cos(theta)));
        beta = 1/(theta*theta) - sin (theta)/(2*theta*(1-cos(theta)));
      }
      matrix3_t rcross; computeCrossMatrix (r, rcross);
      value.segment (0, 3) = alpha*p - .5*rcross*p + beta*r.dot (p)*r;
    }

    template <typename Derived>
      void JlogSE3 (const Transform3f& M,
                    Eigen::MatrixBase<Derived> const& Jlog)
    {
      Eigen::MatrixBase<Derived>& value =
        const_cast<Eigen::MatrixBase<Derived>&> (Jlog);
      const matrix3_t & R = M.rotation();
      const vector3_t & p = M.translation();
      value_type theta;
      vector3_t r;
      logSO3 (R, theta, r);
      matrix3_t Jlog3; JlogSO3 (theta, r, Jlog3);
      value_type alpha, beta, beta_dot_over_theta;
      if (fabs (theta) < 1e-2) {
        alpha = 1 - theta*theta/12 - theta*theta*theta*theta/720;
        beta = 1./12 + theta*theta/720;
        beta_dot_over_theta = 1. / 360.;
      } else {
        alpha = theta*sin(theta)/(2*(1-cos(theta)));
        beta = 1/(theta*theta) - sin (theta)/(2*theta*(1-cos(theta)));
        beta_dot_over_theta = -2/(theta*theta*theta*theta) +
          (theta + sin (theta)) / (2*theta*theta*theta*(1-cos(theta)));
      }
      matrix3_t rcross; computeCrossMatrix (r, rcross);
      matrix3_t V (alpha * matrix3_t::Identity () - .5*rcross +
                   beta * r * r.transpose ());
      value_type rTp (r.dot (p));
      matrix3_t pcross; computeCrossMatrix (p, pcross);
      matrix3_t J ((.5*pcross + (beta_dot_over_theta*rTp)*r*r.transpose ()
                    - (theta*theta*beta_dot_over_theta+2*beta)*p*r.transpose ()
                    + rTp * beta * matrix3_t::Identity ()
                    + beta * r*p.transpose ()) * Jlog3);
      value.block (0, 0, 3, 3) = V * R;
      value.block (0, 3, 3, 3) = J;
      Eigen::Block<Derived> b (value.block (3, 0, 3, 3));
      b.setZero ();
      value.block (3, 3, 3, 3) = Jlog3;
    }

    template <typename Derived1, typename Derived2> void matrixToQuat
    (const Eigen::MatrixBase<Derived1>& M, Eigen::MatrixBase<Derived2> const& q)
    {
      Derived2& _q = const_cast<Derived2&> (q.derived());
      assert (q.size() == 4);
      Eigen::Map<Transform3f::Quaternion> quat (_q.data());
      quat = M;
    }

    template <typename Derived> void se3ToConfig
    (const Transform3f& M, Eigen::MatrixBase<Derived> const& q)
    {
      Derived& _q = const_cast<Derived&> (q.derived());
      assert (q.size() == 7);
      _q.template head<3>() = M.translation();
      matrixToQuat (M.rotation(), _q.template tail<4>());
    }

    /// \}
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_TOOL_HH
