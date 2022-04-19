// Copyright (c) 2017, Joseph Mirabel
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

#ifndef TEST_UTIL_HH
#define TEST_UTIL_HH

#define EIGEN_VECTOR_IS_APPROX(Va, Vb)                    \
  BOOST_CHECK_MESSAGE((Va).isApprox(Vb, test_precision),  \
                      "check " #Va ".isApprox(" #Vb       \
                      ") failed "                         \
                      "[\n"                               \
                          << (Va).transpose() << "\n!=\n" \
                          << (Vb).transpose() << "\n]")

#define SE3CONFIG_IS_APPROX(Va, Vb)                                            \
  {                                                                            \
    BOOST_CHECK_MESSAGE((Va.head<3>()).isApprox(Vb.head<3>(), test_precision), \
                        "check " #Va ".isApprox(" #Vb                          \
                        ") failed "                                            \
                        "[\n"                                                  \
                            << (Va).transpose() << "\n!=\n"                    \
                            << (Vb).transpose() << "\n]");                     \
    if ((Va.tail<4>().array() * Vb.tail<4>().array() > 0.01).any()) {          \
      BOOST_CHECK_MESSAGE(                                                     \
          (Va.tail<4>()).isApprox(Vb.tail<4>(), test_precision),               \
          "check " #Va ".isApprox(" #Vb                                        \
          ") failed "                                                          \
          "[\n"                                                                \
              << (Va).transpose() << "\n!=\n"                                  \
              << (Vb).transpose() << "\n]");                                   \
    } else {                                                                   \
      BOOST_CHECK_MESSAGE(                                                     \
          (Va.tail<4>()).isApprox(-Vb.tail<4>(), test_precision),              \
          "check " #Va ".isApprox(" #Vb                                        \
          ") failed "                                                          \
          "[\n"                                                                \
              << (Va).transpose() << "\n!=\n"                                  \
              << (Vb).transpose() << "\n]");                                   \
    }                                                                          \
  }

#define EIGEN_VECTOR_IS_NOT_APPROX(Va, Vb)                \
  BOOST_CHECK_MESSAGE(!Va.isApprox(Vb, test_precision),   \
                      "check !" #Va ".isApprox(" #Vb      \
                      ") failed "                         \
                      "[\n"                               \
                          << (Va).transpose() << "\n==\n" \
                          << (Vb).transpose() << "\n]")

#define EIGEN_IS_APPROX(matrixA, matrixB)                        \
  BOOST_CHECK_MESSAGE(matrixA.isApprox(matrixB, test_precision), \
                      "check " #matrixA ".isApprox(" #matrixB    \
                      ") failed "                                \
                      "[\n"                                      \
                          << matrixA << "\n!=\n"                 \
                          << matrixB << "\n]")

#define EIGEN_IS_NOT_APPROX(matrixA, matrixB)                     \
  BOOST_CHECK_MESSAGE(!matrixA.isApprox(matrixB, test_precision), \
                      "check !" #matrixA ".isApprox(" #matrixB    \
                      ") failed "                                 \
                      "[\n"                                       \
                          << matrixA << "\n==\n"                  \
                          << matrixB << "\n]")

#define SOLVER_CHECK_SOLVE(expr, expected)                                     \
  {                                                                            \
    hpp::constraints::solver::HierarchicalIterative::Status __status = expr;   \
    BOOST_CHECK_MESSAGE(                                                       \
        __status == hpp::constraints::solver::HierarchicalIterative::expected, \
        "check " #expr " == " #expected " failed ["                            \
            << __status << " != "                                              \
            << hpp::constraints::solver::HierarchicalIterative::expected       \
            << "]");                                                           \
  }

#include <boost/test/included/unit_test.hpp>
#include <hpp/constraints/differentiable-function.hh>
#include <hpp/constraints/fwd.hh>
#include <hpp/constraints/matrix-view.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/extra-config-space.hh>
#include <hpp/pinocchio/joint-collection.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/util.hh>
#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/multibody/model.hpp>

void randomConfig(const hpp::pinocchio::DevicePtr_t& d,
                  hpp::pinocchio::Configuration_t& q) {
  using namespace hpp::constraints;
  size_type extraDim = d->extraConfigSpace().dimension();
  size_type offset = d->configSize() - extraDim;

  q.resize(d->configSize());
  q.head(offset) = ::pinocchio::randomConfiguration(d->model());

  // Shoot extra configuration variables
  for (size_type i = 0; i < extraDim; ++i) {
    value_type lower = d->extraConfigSpace().lower(i);
    value_type upper = d->extraConfigSpace().upper(i);
    value_type range = upper - lower;
    if ((range < 0) || (range == std::numeric_limits<double>::infinity())) {
      std::ostringstream oss;
      oss << "Cannot uniformy sample extra config variable " << i
          << ". min = " << lower << ", max = " << upper << std::endl;
      throw std::runtime_error(oss.str());
    }
    q[offset + i] = lower + (upper - lower) * rand() / RAND_MAX;
  }
}

// This is an ugly fix to make BOOST_CHECK_EQUAL able to print segments_t
// when they are not equal.
namespace std {
std::ostream& operator<<(std::ostream& os,
                         hpp::constraints::BlockIndex::segments_t b) {
  return os << hpp::pretty_print(b);
}
}  // namespace std

class Quadratic : public hpp::constraints::DifferentiableFunction {
 public:
  typedef hpp::shared_ptr<Quadratic> Ptr_t;
  typedef hpp::constraints::DifferentiableFunction DifferentiableFunction;
  typedef hpp::constraints::matrix_t matrix_t;
  typedef hpp::constraints::matrixOut_t matrixOut_t;
  typedef hpp::constraints::vector_t vector_t;
  typedef hpp::constraints::vectorIn_t vectorIn_t;
  typedef hpp::constraints::LiegroupElementRef LiegroupElementRef;
  typedef hpp::constraints::value_type value_type;

  Quadratic(const matrix_t& _A, const vector_t& _b, const value_type& _c)
      : hpp::constraints::DifferentiableFunction(_A.cols(), _A.cols(), 1,
                                                 "Quadratic"),
        A(_A),
        b(_b),
        c(_c) {
    check();
  }

  Quadratic(const matrix_t& _A, const value_type& _c = 0)
      : hpp::constraints::DifferentiableFunction(_A.cols(), _A.cols(), 1,
                                                 "Quadratic"),
        A(_A),
        b(vector_t::Zero(_A.rows())),
        c(_c) {
    check();
  }

  void check() const {
    BOOST_REQUIRE(A.rows() == A.cols());
    BOOST_REQUIRE(A.rows() == b.rows());
  }

  void impl_compute(LiegroupElementRef y, vectorIn_t x) const {
    y.vector()[0] = x.transpose() * A * x + b.dot(x) + c;
  }

  void impl_jacobian(matrixOut_t J, vectorIn_t x) const {
    J.noalias() = 2 * x.transpose() * A + b.transpose();
  }

  matrix_t A;
  vector_t b;
  value_type c;
};

#endif  // TEST_UTIL_HH
