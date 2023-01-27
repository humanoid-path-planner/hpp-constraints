// Copyright (c) 2016, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#ifndef SRC_GENERIC_TRANSFORMATION_HELPER_HH
#define SRC_GENERIC_TRANSFORMATION_HELPER_HH

#include <hpp/constraints/macros.hh>
#include <hpp/constraints/matrix-view.hh>
#include <hpp/constraints/tools.hh>  // for logSO3
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint-collection.hh>
#include <hpp/pinocchio/joint.hh>
#include <pinocchio/multibody/model.hpp>

namespace hpp {
namespace constraints {

namespace {
inline const Transform3f& Id() {
  static const Transform3f id(Transform3f::Identity());
  return id;
}

/// ------- Generic Transform Data ---------------------------------------
template <bool rel>
struct GTOriDataV {};
template <>
struct GTOriDataV<true> {
  Transform3f M;
  value_type theta;
};
template <bool rel>
struct GTOriDataJ {};
template <>
struct GTOriDataJ<true> {
  eigen::matrix3_t Jlog_from1;
};
/// This class contains the data of the GenericTransformation class.
template <bool rel>
struct GTDataBase {
  hpp::pinocchio::DeviceSync device;
  const GenericTransformationModel<rel>& model;
  hpp::pinocchio::DeviceData& ddata() { return device.d(); }

  const JointJacobian_t& J2() { return model.joint2->jacobian(ddata()); }
  const Transform3f& M2() {
    if (model.joint2)
      return model.joint2->currentTransformation(ddata());
    else
      return Id();
  }
  const vector3_t& t2() { return M2().translation(); }
  const matrix3_t& R2() { return M2().rotation(); }

  const JointJacobian_t& J1() {
    return model.getJoint1()->jacobian(this->ddata());
  }
  const Transform3f& M1() {
    return model.getJoint1()->currentTransformation(this->ddata());
  }
  const matrix3_t& R1() { return M1().rotation(); }
  const vector3_t& t1() { return M1().translation(); }

  GTDataBase(const GenericTransformationModel<rel>& m, const DevicePtr_t& d)
      : device(d), model(m) {}
};
template <bool rel, bool pos, bool ori, bool ose3>
struct GTDataV : GTDataBase<rel>, GTOriDataV<ori> {
  enum { ValueSize = (pos ? 3 : 0) + (ori ? (ose3 ? 4 : 3) : 0) };
  typedef Eigen::Matrix<value_type, ValueSize, 1> ValueType;
  ValueType value;

  GTDataV(const GenericTransformationModel<rel>& m, const DevicePtr_t& d)
      : GTDataBase<rel>(m, d) {}
};
/// This class contains the data of the GenericTransformation class.
template <bool rel, bool pos, bool ori, bool ose3>
struct GTDataJ : GTDataV<rel, pos, ori, ose3>, GTOriDataJ<ori> {
  enum {
    ValueSize = GTDataV<rel, pos, ori, ose3>::ValueSize,
    JacobianSize = (pos ? 3 : 0) + (ori ? 3 : 0),
    RowPos = (pos ? 0 : -1),
    RowOri = (ori ? (pos ? 3 : 0) : -1)
  };
  typedef Eigen::Matrix<value_type, JacobianSize, Eigen::Dynamic> JacobianType;
  typedef Eigen::Matrix<value_type, 3, Eigen::Dynamic> matrix3x_t;
  JacobianType jacobian;
  matrix3x_t tmpJac;
  eigen::vector3_t cross1, cross2;

  GTDataJ(const GenericTransformationModel<rel>& m, const DevicePtr_t& d)
      : GTDataV<rel, pos, ori, ose3>(m, d)
  // TODO the two following matrices should be of type Eigen::Map<...>
  // and they should point to some buffer in m.device
  // , jacobian (buffer1, NbRows, m.cols)
  // , tmpJac   (buffer2,      3, m.cols)
  {
    assert(!ose3 || (!ori || m.fullOri));
    if (!m.fullPos || !m.fullOri) jacobian.resize((int)JacobianSize, m.cols);
    cross1.setZero();
    if (m.t2isZero) cross2.setZero();
  }
};

/// ------- Compute log --------------------------------------------------
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

             The expression below assumes the first representation is
             used.
    \param theta angle of rotation \f$R\f$, also \f$\|r\|\f$,
    \param log 3d vector \f$\mathbf{r}\f$,
    \retval Jlog matrix \f$J_{log} (R)\f$.

    \f{align*}
    J_{log} (R) &=&
   \frac{\|\mathbf{r}\|\sin\|\mathbf{r}\|}{2(1-\cos\|\mathbf{r}\|)} I_3 - \frac
   {1}{2}\left[\mathbf{r}\right]_{\times} + (\frac{1}{\|\mathbf{r}\|^2} -
   \frac{\sin\|\mathbf{r}\|}{2\|\mathbf{r}\|(1-\cos\|\mathbf{r}\|)})
   \mathbf{r}\mathbf{r}^T\\
     &=& I_3 -\frac{1}{2}\left[\mathbf{r}\right]_{\times} +
   \left(\frac{2(1-\cos\|\mathbf{r}\|) -
   \|\mathbf{r}\|\sin\|\mathbf{r}\|}{2\|\mathbf{r}\|^2(1-\cos\|\mathbf{r}\|)}\right)\left[\mathbf{r}\right]_{\times}^2
     \f} */
template <typename Derived>
void computeJlog(const value_type& theta, const Eigen::MatrixBase<Derived>& log,
                 matrix3_t& Jlog) {
  if (theta < 1e-6)
    Jlog.setIdentity();
  else {
    // Jlog = alpha I
    const value_type ct = cos(theta), st = sin(theta);
    const value_type st_1mct = st / (1 - ct);

    Jlog.setZero();
    Jlog.diagonal().setConstant(theta * st_1mct);

    // Jlog += -r_{\times}/2
    Jlog(0, 1) = log(2);
    Jlog(1, 0) = -log(2);
    Jlog(0, 2) = -log(1);
    Jlog(2, 0) = log(1);
    Jlog(1, 2) = log(0);
    Jlog(2, 1) = -log(0);
    Jlog /= 2;

    const value_type alpha = 1 / (theta * theta) - st_1mct / (2 * theta);
    Jlog.noalias() += alpha * log * log.transpose();
  }
}

typedef JointJacobian_t::ConstNRowsBlockXpr<3>::Type HalfJacobian_t;
inline HalfJacobian_t omega(const JointJacobian_t& j) {
  return j.bottomRows<3>();
}
inline HalfJacobian_t trans(const JointJacobian_t& j) { return j.topRows<3>(); }

static inline size_type size(std::vector<bool> mask) {
  size_type res = 0;
  for (std::vector<bool>::iterator it = mask.begin(); it != mask.end(); ++it) {
    if (*it) ++res;
  }
  return res;
}

template <bool flag /* false */>
struct unary {
  template <bool rel, bool pos, bool ose3>
  static inline void log(GTDataV<rel, pos, flag, ose3>&) {}
  template <bool rel, bool pos, bool ose3>
  static inline void Jlog(GTDataJ<rel, pos, flag, ose3>&) {}
};
template <>
struct unary<true> {
  template <bool rel, bool pos, bool ose3>
  static inline void log(GTDataV<rel, pos, true, ose3>& d) {
    if (ose3) {
      matrixToQuat(d.M.rotation(), d.value.template tail<4>());
    } else {
      logSO3(d.M.rotation(), d.theta, d.value.template tail<3>());
      hppDnum(info, "theta=" << d.theta);
    }
  }
  template <bool rel, bool pos, bool ose3>
  static inline void Jlog(GTDataJ<rel, pos, true, ose3>& d) {
    if (ose3) {
      d.Jlog_from1 = d.model.F2inJ2.rotation().transpose() * d.R2().transpose();
      if (rel && d.model.getJoint1()) d.Jlog_from1 *= d.R1();
    } else {
      computeJlog(d.theta, d.value.template tail<3>(), d.Jlog_from1);
      hppDnum(info, "Jlog_: " << d.Jlog_from1);
      if (!d.model.R1isID)
        d.Jlog_from1 *= d.model.F1inJ1.rotation().transpose();
    }
  }
};

template <bool ori, typename Data, typename Derived>
void assign_if(bool cond, Data& d, matrixOut_t J,
               const Eigen::MatrixBase<Derived>& rhs,
               const size_type& startRow) {
  const int& rowCache = (ori ? Data::RowOri : Data::RowPos);
  if (cond)
    d.jacobian.template middleRows<3>(rowCache).noalias() = rhs;
  else
    J.template middleRows<3>(startRow).leftCols(d.model.cols).noalias() = rhs;
}

/// ------- Compute jacobian ---------------------------------------------
template <bool lflag /*rel*/, bool rflag /*false*/>
struct binary {
  // the first template allow us to consider relative transformation as
  // absolute when joint1 is NULL, at run time
  template <bool rel, bool pos, bool ose3>
  static inline void Jorientation(GTDataJ<rel, pos, rflag, ose3>&,
                                  matrixOut_t) {}
  template <bool rel, bool ori, bool ose3>
  static inline void Jtranslation(GTDataJ<rel, rflag, ori, ose3>&,
                                  matrixOut_t) {}
};
template <>
struct binary<false, true>  // Absolute
{
  template <bool rel, bool pos, bool ose3>
  static inline void Jorientation(GTDataJ<rel, pos, true, ose3>& d,
                                  matrixOut_t J) {
    assert(!ose3 || d.model.fullOri);
    assign_if<true>(!(ose3 || d.model.fullOri), d, J,
                    (d.Jlog_from1 * d.R2()) * omega(d.J2()), d.model.rowOri);
  }
  template <bool rel, bool ori, bool ose3>
  static inline void Jtranslation(GTDataJ<rel, true, ori, ose3>& d,
                                  matrixOut_t J) {
    const JointJacobian_t& J2(d.J2());
    const matrix3_t& R2(d.R2());
    const matrix3_t& R1inJ1(d.model.F1inJ1.rotation());

    // hpp-model: J = 1RT* ( 0Jt2 - [ 0R2 2t* ]x 0Jw2 )
    // pinocchio: J = 1RT* ( 0R2 2Jt2 - [ 0R2 2t* ]x 0R2 2Jw2 )
    if (!d.model.t2isZero) {
      d.tmpJac.noalias() = (R2.colwise().cross(d.cross2)) * omega(J2);
      d.tmpJac.noalias() += R2 * trans(J2);
      if (d.model.R1isID) {
        assign_if<false>(!d.model.fullPos, d, J, d.tmpJac, 0);
      } else {  // Generic case
        assign_if<false>(!d.model.fullPos, d, J, R1inJ1.transpose() * d.tmpJac,
                         0);
      }
    } else {
      if (d.model.R1isID)
        assign_if<false>(!d.model.fullPos, d, J, R2 * trans(J2), 0);
      else
        assign_if<false>(!d.model.fullPos, d, J,
                         (R1inJ1.transpose() * R2) * trans(J2), 0);
    }
  }
};
template <>
struct binary<true, true>  // Relative
{
  template <bool pos, bool ose3>
  static inline void Jorientation(GTDataJ<true, pos, true, ose3>& d,
                                  matrixOut_t J) {
    d.tmpJac.noalias() = -omega(d.J1());
    if (d.model.joint2)
      d.tmpJac.noalias() += (d.R1().transpose() * d.R2()) * omega(d.J2());
    assert(!ose3 || d.model.fullOri);
    assign_if<true>(!(ose3 || d.model.fullOri), d, J, d.Jlog_from1 * d.tmpJac,
                    d.model.rowOri);
  }
  template <bool ori, bool ose3>
  static inline void Jtranslation(GTDataJ<true, true, ori, ose3>& d,
                                  matrixOut_t J) {
    const JointJacobian_t& J1(d.J1());
    const matrix3_t& R1(d.R1());
    const matrix3_t& R2(d.R2());
    const matrix3_t& R1inJ1(d.model.F1inJ1.rotation());

    // J = 1RT* 0RT1 ( A + B )
    // hpp-model:
    // A = [ 0t2 - 0t1 0R2 2t* ]x 0Jw1
    // B = ( 0Jt2 - 0Jt1 - [ 0R2 2t* ]x 0Jw2 )
    // pinocchio:
    // A = [ 0t2 - 0t1 0R2 2t* ]x 0R1 1Jw1
    // B = ( 0R2 2Jt2 - 0R1 1Jt1 - [ 0R2 2t* ]x 0R2 2Jw2 )
    d.tmpJac.noalias() =
        (-R1.transpose() * R1.colwise().cross(d.cross1)) * omega(J1);  // A
    if (d.model.joint2)
      d.tmpJac.noalias() += (R1.transpose() * R2) * trans(d.J2());  // B1
    d.tmpJac.noalias() -= trans(J1);                                // B2
    if (!d.model.t2isZero && d.model.joint2)
      d.tmpJac.noalias() +=
          R1.transpose() * R2.colwise().cross(d.cross2) * omega(d.J2());  // B3
    if (d.model.R1isID)
      assign_if<false>(!d.model.fullPos, d, J, d.tmpJac, 0);
    else
      assign_if<false>(!d.model.fullPos, d, J, R1inJ1.transpose() * d.tmpJac,
                       0);
  }
};

/// ------- Compute relative transform -----------------------------------
template <bool compileTimeRel /* false */, bool ori /* false */>
struct relativeTransform {
  template <bool runtimeRel, bool pos, bool ose3>
  static inline void run(GTDataV<runtimeRel, pos, false, ose3>& d) {
    using hpp::pinocchio::LiegroupElement;
    using hpp::pinocchio::LiegroupSpace;
    // There is no joint1
    const Transform3f& M2 = d.M2();
    d.value.noalias() = M2.act(d.model.F2inJ2.translation());
#ifndef NDEBUG
    if (ose3) {
      hpp::pinocchio::vector4_t quat;
      size_type n = d.value.size();
      assert(n >= 4);
      quat << d.value[n - 4], d.value[n - 3], d.value[n - 2], d.value[n - 1];
      assert(hpp::pinocchio::checkNormalized(
          LiegroupElement(quat, LiegroupSpace::SO3())));
    }
#endif
    if (!d.model.t1isZero) d.value.noalias() -= d.model.F1inJ1.translation();
#ifndef NDEBUG
    if (ose3) {
      hpp::pinocchio::vector4_t quat;
      size_type n = d.value.size();
      assert(n >= 4);
      quat << d.value[n - 4], d.value[n - 3], d.value[n - 2], d.value[n - 1];
      assert(hpp::pinocchio::checkNormalized(
          LiegroupElement(quat, LiegroupSpace::SO3())));
    }
#endif
    if (!d.model.R1isID)
      d.value.applyOnTheLeft(d.model.F1inJ1.rotation().transpose());
#ifndef NDEBUG
    if (ose3) {
      hpp::pinocchio::vector4_t quat;
      size_type n = d.value.size();
      assert(n >= 4);
      quat << d.value[n - 4], d.value[n - 3], d.value[n - 2], d.value[n - 1];
      assert(hpp::pinocchio::checkNormalized(
          LiegroupElement(quat, LiegroupSpace::SO3())));
    }
#endif
  }
};
template <>
struct relativeTransform<false, true> {
  template <bool runtimeRel, bool pos, bool ose3>
  static inline void run(GTDataV<runtimeRel, pos, true, ose3>& d) {
    const Transform3f& M2 = d.M2();
    d.M = d.model.F1inJ1.actInv(M2 * d.model.F2inJ2);
    if (pos) d.value.template head<3>().noalias() = d.M.translation();
  }
};
template <>
struct relativeTransform<true, true> {
  template <bool pos, bool ose3>
  static inline void run(GTDataV<true, pos, true, ose3>& d) {
    if (d.model.joint1 == NULL) {
      // runtime absolute reference.
      relativeTransform<false, true>::run(d);
      return;
    }
    const Transform3f& M1 = d.M1();
    const Transform3f& M2 = d.M2();
    d.M = d.model.F1inJ1.actInv(M1.actInv(M2 * d.model.F2inJ2));
    if (pos) d.value.template head<3>().noalias() = d.M.translation();
  }
};
template <>
struct relativeTransform<true, false> {
  template <bool pos, bool ose3>
  static inline void run(GTDataV<true, pos, false, ose3>& d) {
    if (d.model.joint1 == NULL) {
      // runtime absolute reference.
      relativeTransform<false, false>::run(d);
      return;
    }
    const Transform3f& M2 = d.M2();
    const Transform3f& M1 = d.M1();
    d.value.noalias() = M2.act(d.model.F2inJ2.translation()) - M1.translation();
    d.value.applyOnTheLeft(M1.rotation().transpose());

    if (!d.model.t1isZero) d.value.noalias() -= d.model.F1inJ1.translation();
    if (!d.model.R1isID)
      d.value.applyOnTheLeft(d.model.F1inJ1.rotation().transpose());
  }
};

template <bool rel, bool pos, bool ori, bool ose3>
struct compute {
  static inline void error(GTDataV<rel, pos, ori, ose3>& d) {
    using hpp::pinocchio::LiegroupElement;
    using hpp::pinocchio::LiegroupSpace;
    relativeTransform<rel, ori>::run(d);
    unary<ori>::log(d);
#ifndef NDEBUG
    if (ose3) {
      hpp::pinocchio::vector4_t quat;
      size_type n = d.value.size();
      assert(n >= 4);
      quat << d.value[n - 4], d.value[n - 3], d.value[n - 2], d.value[n - 1];
      assert(hpp::pinocchio::checkNormalized(
          LiegroupElement(quat, LiegroupSpace::SO3())));
    }
#endif
  }

  static inline void jacobian(GTDataJ<rel, pos, ori, ose3>& d,
                              matrixOut_t jacobian,
                              const std::vector<bool>& mask) {
    const Transform3f& M2 = d.M2();
    const vector3_t& t2inJ2(d.model.F2inJ2.translation());
    const vector3_t& t2(M2.translation());
    const matrix3_t& R2(M2.rotation());

    if (!d.model.t2isZero) d.cross2.noalias() = R2 * t2inJ2;

    unary<ori>::Jlog(d);

    // rel:           relative known at compile time
    // d.getJoint1(): relative known at run time
    if (rel && d.model.getJoint1()) {
      const Transform3f& M1 = d.M1();
      const vector3_t& t1(M1.translation());
      d.cross1.noalias() = d.cross2 + t2 - t1;
      binary<rel, pos>::Jtranslation(d, jacobian);
      binary<rel, ori>::Jorientation(d, jacobian);
    } else {
      d.cross1.noalias() = d.cross2 + t2;
      binary<false, pos>::Jtranslation(d, jacobian);
      binary<false, ori>::Jorientation(d, jacobian);
    }

    // Copy necessary rows.
    size_type index = 0;
    const size_type lPos = (pos ? 3 : 0), lOri = (ori ? 3 : 0);
    if (!d.model.fullPos) {
      for (size_type i = 0; i < lPos; ++i) {
        if (mask[i]) {
          jacobian.row(index).leftCols(d.model.cols).noalias() =
              d.jacobian.row(i);
          ++index;
        }
      }
    } else
      index = lPos;
    if (!d.model.fullOri) {
      for (size_type i = lPos; i < lPos + lOri; ++i) {
        if (mask[i]) {
          jacobian.row(index).leftCols(d.model.cols).noalias() =
              d.jacobian.row(i);
          ++index;
        }
      }
    }
    jacobian.rightCols(jacobian.cols() - d.model.cols).setZero();
  }
};

void setActiveParameters(const DevicePtr_t& robot, const JointConstPtr_t& j1,
                         const JointConstPtr_t& j2, ArrayXb& activeParameters,
                         ArrayXb& activeDerivativeParameters) {
  typedef ::pinocchio::JointIndex JointIndex;

  const pinocchio::Model& model = robot->model();
  const JointIndex id1 = (j1 ? j1->index() : 0), id2 = (j2 ? j2->index() : 0);
  JointIndex i1 = id1, i2 = id2;

  std::vector<JointIndex> from1, from2;
  while (i1 != i2) {
    JointIndex i;
    if (i1 > i2) {
      i = i1;
      i1 = model.parents[i1];
    } else /* if (i1 < i2) */ {
      i = i2;
      i2 = model.parents[i2];
    }
    if (i > 0) {
      activeParameters.segment(model.joints[i].idx_q(), model.joints[i].nq())
          .setConstant(true);
      activeDerivativeParameters
          .segment(model.joints[i].idx_v(), model.joints[i].nv())
          .setConstant(true);
    }
  }
  assert(i1 == i2);
}
}  // namespace
}  // namespace constraints
}  // namespace hpp

#endif  // SRC_GENERIC_TRANSFORMATION_HELPER_HH
