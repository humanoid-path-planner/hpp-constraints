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

#include <hpp/constraints/fwd.hh>

namespace hpp {
  namespace constraints {
    template <typename Derived>
      inline void computeLog (const matrix3_t& Rerror,
          value_type& theta, Eigen::MatrixBase<Derived> const& result)
    {
      Eigen::MatrixBase<Derived>& value = const_cast<Eigen::MatrixBase<Derived>&> (result);
      const value_type PI = ::boost::math::constants::pi<value_type>();
      const value_type tr = Rerror.trace();
      if (tr > 3)       theta = 0; // acos((3-1)/2)
      else if (tr < -1) theta = PI; // acos((-1-1)/2)
      else              theta = acos ((tr - 1)/2);
      assert (theta == theta);
      // From runs of tests/logarithm.cc: 1e-6 is too small.
      if (theta < PI - 1e-2) {
        const value_type t = ((theta > 1e-6)? theta / sin(theta) : 1) / 2;
        value(0) = t * (Rerror (2, 1) - Rerror (1, 2));
        value(1) = t * (Rerror (0, 2) - Rerror (2, 0));
        value(2) = t * (Rerror (1, 0) - Rerror (0, 1));
      } else {
        // 1e-2: A low value is not required since the computation is
        // using explicit formula. However, the precision of this method
        // is the square root of the precision with the antisymmetric
        // method (Nominal case).
        const value_type cphi = cos(theta - PI);
        const value_type beta  = theta*theta / ( 1 + cphi );
        const value_type tmp0 = (Rerror (0, 0) + cphi) * beta;
        const value_type tmp1 = (Rerror (1, 1) + cphi) * beta;
        const value_type tmp2 = (Rerror (2, 2) + cphi) * beta;
        value(0) = (Rerror (2, 1) > Rerror (1, 2) ? 1 : -1 ) * (tmp0 > 0 ? sqrt(tmp0) : 0);
        value(1) = (Rerror (0, 2) > Rerror (2, 0) ? 1 : -1 ) * (tmp1 > 0 ? sqrt(tmp1) : 0);
        value(2) = (Rerror (1, 0) > Rerror (0, 1) ? 1 : -1 ) * (tmp2 > 0 ? sqrt(tmp2) : 0);
      }
    }


    /// Compute jacobian of function log of rotation matrix in SO(3)
    ///
    /// Let us consider a matrix
    /// \f$R=\exp \left[\mathbf{r}\right]_{\times}\in SO(3)\f$.
    /// This functions computes the Jacobian of the function from
    /// \f$SO(3)\f$ into \f$\mathbf{R}^3\f$ that maps \f$R\f$ to
    /// \f$\mathbf{r}\f$. In other words,
    /// \f{equation*}
    /// \dot{\mathbf{r}} = J_{log}(R)\ \omega\,\,\,\mbox{with}\,\,\,
    /// \dot {R} = \left[\omega\right]_{\times} R
    /// \f}
    /// \warning Two representations of the angular velocity \f$\omega\f$ are
    ///          possible:
    ///          \li \f$\dot{R} = \left[\omega\right]_{\times}R\f$ or
    ///          \li \f$\dot{R} = R\left[\omega\right]_{\times}\f$.
    ///
    ///          The expression below is different with the second
    ///          representation.
    /// \param theta angle of rotation \f$R\f$, also \f$\|r\|\f$,
    /// \param log 3d vector \f$\mathbf{r}\f$,
    /// \retval Jlog matrix \f$J_{log} (R)\f$.
    ///
    /// \f{align*}
    /// J_{log} (R) &=& \frac{\|\mathbf{r}\|\sin\|\mathbf{r}\|}{2(1-\cos\|\mathbf{r}\|)} I_3 - \frac {1}{2}\left[\mathbf{r}\right]_{\times} + (\frac{1}{\|\mathbf{r}\|^2} - \frac{\sin\|\mathbf{r}\|}{2\|\mathbf{r}\|(1-\cos\|\mathbf{r}\|)}) \mathbf{r}\mathbf{r}^T\\
    ///  &=& I_3 -\frac{1}{2}\left[\mathbf{r}\right]_{\times} +  \left(\frac{2(1-\cos\|\mathbf{r}\|) - \|\mathbf{r}\|\sin\|\mathbf{r}\|}{2\|\mathbf{r}\|^2(1-\cos\|\mathbf{r}\|)}\right)\left[\mathbf{r}\right]_{\times}^2
    /// \f}
    template <typename Derived>
      void computeJlog (const value_type& theta, const Eigen::MatrixBase<Derived>& log, matrix3_t& Jlog)
    {
      if (theta < 1e-6)
        Jlog.setIdentity();
      else {
        // Jlog = alpha I
        const value_type ct = cos(theta), st = sin(theta);
        const value_type st_1mct = st/(1-ct);

        Jlog.setZero ();
        Jlog.diagonal().setConstant (theta*st_1mct);

        // Jlog += -r_{\times}/2
        Jlog(0,1) =  log(2); Jlog(1,0) = -log(2);
        Jlog(0,2) = -log(1); Jlog(2,0) =  log(1);
        Jlog(1,2) =  log(0); Jlog(2,1) = -log(0);
        Jlog /= 2;

        const value_type alpha = 1/(theta*theta) - st_1mct/(2*theta);
        Jlog.noalias() += alpha * log * log.transpose ();
      }
    }

    template < typename VectorType, typename MatrixType >
    static void computeCrossMatrix (const VectorType& v, MatrixType& m)
    {
      m.diagonal ().setZero ();
      m (0,1) = -v [2]; m (1,0) = v [2];
      m (0,2) = v [1]; m (2,0) = -v [1];
      m (1,2) = -v [0]; m (2,1) = v [0];
    }
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_TOOL_HH
