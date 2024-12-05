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

#define BOOST_TEST_MODULE hpp_constraints
#include <boost/test/included/unit_test.hpp>
#include <hpp/constraints/tools.hh>
#include <pinocchio/fwd.hpp>
#include <pinocchio/spatial/explog.hpp>
#include <pinocchio/spatial/motion.hpp>

using namespace hpp::constraints;

matrix3_t exponential(const vector3_t& aa) {
  matrix3_t R, xCross;
  xCross.setZero();
  xCross(1, 0) = +aa(2);
  xCross(0, 1) = -aa(2);
  xCross(2, 0) = -aa(1);
  xCross(0, 2) = +aa(1);
  xCross(2, 1) = +aa(0);
  xCross(1, 2) = -aa(0);
  R.setIdentity();
  value_type theta = aa.norm();
  if (theta < 1e-6) {
    R += xCross;
    R += 0.5 * xCross.transpose() * xCross;
  } else {
    R += sin(theta) / theta * xCross;
    R += 2 * std::pow(sin(theta / 2), 2) / std::pow(theta, 2) * xCross * xCross;
  }
  return R;
}

vector3_t log(const matrix3_t R) {
  vector3_t log;
  value_type theta;
  logSO3(R, theta, log);
  return log;
}

bool check(const vector3_t& aa, const value_type eps = -1) {
  vector3_t res = log(exponential(aa));
  bool ret;
  if (eps < 0)
    ret = aa.isApprox(res);
  else
    ret = aa.isApprox(res, eps);
  if (!ret) {
    std::cout << "x = " << aa << std::endl;
    std::cout << "exp(x) = " << exponential(aa) << std::endl;
    std::cout << "log(exp(x)) = " << res << std::endl;
    std::cout << "(x - log(exp(x))) / norm(x) = "
              << (aa - res).derived().transpose() / aa.norm() << std::endl;
    std::cout << "eps = "
              << (eps < 0 ? Eigen::NumTraits<value_type>::dummy_precision()
                          : eps)
              << std::endl;
  }
  return ret;
}

BOOST_AUTO_TEST_CASE(logarithm) {
  // Width of the domain considered as close to rotation of Pi. It should be the
  // same as the value in the logarithm computation.
  // When in this domain, the precision of the computation is done on the square
  // of the angle axis. So we only have precision eps defined as follow.
  const value_type dlUB = 1e-2;
  const value_type eps = sqrt(Eigen::NumTraits<value_type>::epsilon());

  BOOST_CHECK(check(vector3_t(0, 0, 0)));
  BOOST_CHECK(check(vector3_t(1, 0, 0)));
  BOOST_CHECK(check(vector3_t(0, 1, 0)));
  BOOST_CHECK(check(vector3_t(0, 0, 1)));
  BOOST_CHECK(check(vector3_t(1, 1, 0)));
  BOOST_CHECK(check(vector3_t(0, 1, 1)));
  BOOST_CHECK(check(vector3_t(1, 0, 1)));
  BOOST_CHECK(check(vector3_t(1, 1, 1)));

  BOOST_CHECK(check(M_PI * vector3_t(1, 0, 0)));
  BOOST_CHECK(check(M_PI * vector3_t(0, 1, 0)));
  BOOST_CHECK(check(M_PI * vector3_t(0, 0, 1)));
  BOOST_CHECK(check(M_PI / sqrt(2) * vector3_t(1, 1, 0)));
  BOOST_CHECK(check(M_PI / sqrt(2) * vector3_t(0, 1, 1)));
  BOOST_CHECK(check(M_PI / sqrt(2) * vector3_t(1, 0, 1)));

  BOOST_CHECK(check(M_PI / sqrt(3) * vector3_t(1, -1, 1), eps));
  BOOST_CHECK(check(M_PI / sqrt(3) * vector3_t(1, 1, 1), eps));

  BOOST_CHECK(check(1e-7 * vector3_t(1, 1, 1)));
  BOOST_CHECK(check((M_PI - 0.1 * dlUB) / sqrt(3) * vector3_t(1, -1, 1), eps));
  BOOST_CHECK(check((M_PI - 0.1 * dlUB) / sqrt(3) * vector3_t(1, 1, 1), eps));

  BOOST_CHECK(check((M_PI - 1.001 * dlUB) / sqrt(3) * vector3_t(1, 1, 1), eps));
  BOOST_CHECK(check((M_PI - 0.999 * dlUB) / sqrt(3) * vector3_t(1, 1, 1)));
}

BOOST_AUTO_TEST_CASE(Jlog_SO3) {
  value_type lower(-1), upper(1);
  for (size_type i = 0; i < 100; ++i) {
    // Generate random rotation matrix
    vector3_t r0;
    r0[0] = lower + (upper - lower) * rand() / RAND_MAX;
    r0[1] = lower + (upper - lower) * rand() / RAND_MAX;
    r0[2] = lower + (upper - lower) * rand() / RAND_MAX;
    r0.normalize();
    r0 *= 3.14 * rand() / RAND_MAX;
    if (i == 0) r0.setZero();
    matrix3_t R0(exponential(r0));
    vector3_t omega;
    // \dot{R} = R0 [\omega]_{\times}
    // R (dt) = R0 \exp (dt [\omega]_{\times})
    value_type dt(1e-6);
    for (size_type i = 0; i < 3; ++i) {
      omega.setZero();
      omega[i] = 1;
      matrix3_t R(R0 * exponential(dt * omega));
      value_type theta;
      vector3_t r;
      logSO3(R, theta, r);
      matrix3_t Jlog;
      JlogSO3(theta, r, Jlog);
      BOOST_CHECK(((r - r0) / dt - Jlog * omega).squaredNorm() < 1e-8);
    }
  }
}

BOOST_AUTO_TEST_CASE(Jlog_SE3) {
  for (size_type i = 0; i < 100; ++i) {
    // Generate random rigid body motion
    Transform3s M0;
    M0.setRandom();
    if (i == 0) M0.setIdentity();
    vector6_t log0, log;
    value_type dt(1e-6);
    logSE3(M0, log0);
    matrix6_t Jlog;
    JlogSE3(M0, Jlog);
    for (size_type i = 0; i < 6; ++i) {
      vector6_t v;
      v.setZero();
      v[i] = 1;
      ::pinocchio::Motion nu(v);
      Transform3s M(M0 * ::pinocchio::exp6(dt * nu));
      logSE3(M, log);
      std::cout << "(log-log0)/dt=" << ((log - log0) / dt).transpose()
                << std::endl;
      std::cout << "Jlog * v=     " << (Jlog * v).transpose() << std::endl;
      BOOST_CHECK(((log - log0) / dt - Jlog * v).squaredNorm() < 1e-8);
    }
  }
}
