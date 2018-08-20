// Copyright (c) 2017, Joseph Mirabel
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

#ifndef TEST_UTIL_HH
# define TEST_UTIL_HH

#define EIGEN_VECTOR_IS_APPROX(Va, Vb)                                         \
  BOOST_CHECK_MESSAGE((Va).isApprox(Vb, test_precision),                       \
      "check " #Va ".isApprox(" #Vb ") failed "                                \
      "[\n" << (Va).transpose() << "\n!=\n" << (Vb).transpose() << "\n]")

#define EIGEN_VECTOR_IS_NOT_APPROX(Va, Vb)                                     \
  BOOST_CHECK_MESSAGE(!Va.isApprox(Vb, test_precision),                        \
      "check !" #Va ".isApprox(" #Vb ") failed "                               \
      "[\n" << (Va).transpose() << "\n==\n" << (Vb).transpose() << "\n]")

#define EIGEN_IS_APPROX(matrixA, matrixB)                                      \
  BOOST_CHECK_MESSAGE(matrixA.isApprox(matrixB, test_precision),               \
      "check " #matrixA ".isApprox(" #matrixB ") failed "                      \
      "[\n" << matrixA << "\n!=\n" << matrixB << "\n]")

#define EIGEN_IS_NOT_APPROX(matrixA, matrixB)                                  \
  BOOST_CHECK_MESSAGE(!matrixA.isApprox(matrixB, test_precision),              \
      "check !" #matrixA ".isApprox(" #matrixB ") failed "                     \
      "[\n" << matrixA << "\n==\n" << matrixB << "\n]")

#define SOLVER_CHECK_SOLVE(expr,expected)                                      \
  {                                                                            \
    hpp::constraints::solver::HierarchicalIterative::Status __status = expr; \
    BOOST_CHECK_MESSAGE(__status == hpp::constraints::solver::HierarchicalIterative::expected, \
        "check " #expr " == " #expected " failed ["                            \
                        << __status << " != " << hpp::constraints::solver::HierarchicalIterative::expected << "]"); \
  }


#include <pinocchio/multibody/model.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/extra-config-space.hh>
#include <hpp/pinocchio/joint.hh>

#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/differentiable-function.hh>

bool saturate (const hpp::pinocchio::DevicePtr_t& robot,
    hpp::pinocchio::vectorIn_t q, hpp::pinocchio::vectorOut_t qSat,
    Eigen::VectorXi& sat)
{
  using hpp::pinocchio::size_type;
  bool ret = false;
  qSat = q;
  const se3::Model& model = robot->model();

  for (std::size_t i = 1; i < model.joints.size(); ++i) {
    const size_type nq = model.joints[i].nq();
    const size_type nv = model.joints[i].nv();
    const size_type idx_q = model.joints[i].idx_q();
    const size_type idx_v = model.joints[i].idx_v();
    for (size_type j = 0; j < nq; ++j) {
      const size_type iq = idx_q + j;
      const size_type iv = idx_v + std::min(j,nv-1);
      if        (q[iq] >= model.upperPositionLimit[iq]) {
        sat[iv] =  1;
        qSat [iq] = model.upperPositionLimit[iq];
        ret = true;
      } else if (q[iq] <= model.lowerPositionLimit[iq]) {
        sat[iv] = -1;
        qSat [iq] = model.lowerPositionLimit[iq];
        ret = true;
      } else
        sat[iv] =  0;
    }
  }

  const hpp::pinocchio::ExtraConfigSpace& ecs = robot->extraConfigSpace();
  const size_type& d = ecs.dimension();

  for (size_type k = 0; k < d; ++k) {
    const size_type iq = model.nq + k;
    const size_type iv = model.nv + k;
    if        (q[iq] >= ecs.upper(k)) {
      sat[iv] =  1;
      ret = true;
    } else if (q[iq] <= ecs.lower(k)) {
      sat[iv] = -1;
      ret = true;
    } else
      sat[iv] =  0;
  }
  return ret;
}

template <int Lower, int Upper>
bool simpleSaturation (hpp::constraints::vectorIn_t x,
                       hpp::constraints::vectorOut_t xSat,
                       Eigen::VectorXi& sat)
{
  bool ret = false;
  xSat = x;
  for (hpp::constraints::size_type i = 0; i < x.size(); ++i) {
    if (x[i] <= Lower) {
      sat[i] = -1;
      xSat [i] = Lower;
      ret = true;
    }
    else if (x[i] >= Upper) {
      sat[i] =  1;
      xSat [i] = Upper;
      ret = true;
    } else {
      sat[i] = 0;
    }
  }
  return ret;
}

// This is an ugly fix to make BOOST_CHECK_EQUAL able to print segments_t
// when they are not equal.
namespace std {
  std::ostream& operator<< (std::ostream& os, hpp::constraints::BlockIndex::segments_t b)
  {
    return os << hpp::pretty_print (b);
  }
}

class Quadratic : public hpp::constraints::DifferentiableFunction
{
  public:
    typedef boost::shared_ptr<Quadratic> Ptr_t;
    typedef hpp::constraints::DifferentiableFunction DifferentiableFunction;
    typedef hpp::constraints::matrix_t matrix_t;
    typedef hpp::constraints::matrixOut_t matrixOut_t;
    typedef hpp::constraints::vector_t vector_t;
    typedef hpp::constraints::vectorIn_t vectorIn_t;
    typedef hpp::constraints::LiegroupElement LiegroupElement;
    typedef hpp::constraints::value_type value_type;

    Quadratic (const matrix_t& _A, const vector_t& _b, const value_type& _c)
      : hpp::constraints::DifferentiableFunction (_A.cols(), _A.cols(), 1, "Quadratic"),
      A (_A), b (_b), c(_c)
    {
      check();
    }

    Quadratic (const matrix_t& _A, const value_type& _c = 0)
      : hpp::constraints::DifferentiableFunction (_A.cols(), _A.cols(), 1, "Quadratic"),
      A (_A), b (vector_t::Zero(_A.rows())), c(_c)
    {
      check();
    }

    void check () const
    {
      BOOST_REQUIRE (A.rows() == A.cols());
      BOOST_REQUIRE (A.rows() == b.rows());
    }

    void impl_compute (LiegroupElement& y, vectorIn_t x) const
    {
      y.vector()[0] = x.transpose() * A * x + b.dot(x) + c;
    }

    void impl_jacobian (matrixOut_t J, vectorIn_t x) const
    {
      J.noalias() = 2 * x.transpose() * A + b.transpose();
    }

    matrix_t A;
    vector_t b;
    value_type c;
};

#endif // TEST_UTIL_HH
