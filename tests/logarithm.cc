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

// #include <hpp/model/device.hh>
// #include <hpp/model/joint.hh>
// #include <hpp/model/configuration.hh>
// #include <hpp/model/object-factory.hh>

// #include "hpp/constraints/position.hh"
// #include "hpp/constraints/orientation.hh"
// #include "hpp/constraints/transformation.hh"
// #include "hpp/constraints/relative-position.hh"
// #include "hpp/constraints/relative-orientation.hh"
// #include "hpp/constraints/relative-transformation.hh"
// #include "hpp/constraints/generic-transformation.hh"
// #include "hpp/constraints/tools.hh"

#define BOOST_TEST_MODULE hpp_constraints
#include <boost/test/included/unit_test.hpp>

#include <../src/generic-transformation.cc>

using namespace hpp::constraints;

matrix3_t exponential (const vector3_t& aa)
{
  matrix3_t R, xCross;
  xCross.setZero();
  xCross(1, 0) = + aa(2); xCross(0, 1) = - aa(2);
  xCross(2, 0) = - aa(1); xCross(0, 2) = + aa(1);
  xCross(2, 1) = + aa(0); xCross(1, 2) = - aa(0);
  R.setIdentity();
  value_type theta = aa.norm();
  if (theta < 1e-6) {
    R += xCross;
    R += 0.5 * xCross.transpose() * xCross;
  } else {
    R += sin(theta) / theta * xCross;
    R += 2 * std::pow(sin(theta/2),2) / std::pow(theta,2) * xCross.transpose() * xCross;
  }
  return R;
}

vector3_t log (const matrix3_t R)
{
  GenericTransformationData<false, false, true> data(0);
  data.M = R;
  unary<true>::log(data);
  return data.value;
}

bool check (const vector3_t& aa, const value_type eps = -1)
{
  vector3_t res = log (exponential(aa));
  bool ret;
  if (eps < 0) ret = aa.isApprox(res);
  else         ret = aa.isApprox(res, eps);
  if (!ret) {
    std::cout << "x = " << aa << std::endl;
    std::cout << "exp(x) = " << exponential(aa) << std::endl;
    std::cout << "log(exp(x)) = " << res << std::endl;
    std::cout << "(x - log(exp(x))) / norm(x) = " << (aa-res).derived().transpose()/aa.norm() << std::endl;
    std::cout << "eps = " << (eps < 0 ? Eigen::NumTraits<value_type>::dummy_precision() : eps) << std::endl;
  }
  return ret;
}

BOOST_AUTO_TEST_CASE (logarithm) {
  // Width of the domain considered as close to rotation of Pi. It should be the
  // same as the value in the logarithm computation.
  // When in this domain, the precision of the computation is done on the square
  // of the angle axis. So we only have precision eps defined as follow.
  const value_type dlUB = 1e-2;
  const value_type eps = sqrt(Eigen::NumTraits<value_type>::epsilon());

  BOOST_CHECK(check (vector3_t (0,0,0)));
  BOOST_CHECK(check (vector3_t (1,0,0)));
  BOOST_CHECK(check (vector3_t (0,1,0)));
  BOOST_CHECK(check (vector3_t (0,0,1)));
  BOOST_CHECK(check (vector3_t (1,1,0)));
  BOOST_CHECK(check (vector3_t (0,1,1)));
  BOOST_CHECK(check (vector3_t (1,0,1)));
  BOOST_CHECK(check (vector3_t (1,1,1)));

  BOOST_CHECK(check (M_PI * vector3_t (1,0,0)));
  BOOST_CHECK(check (M_PI * vector3_t (0,1,0)));
  BOOST_CHECK(check (M_PI * vector3_t (0,0,1)));
  BOOST_CHECK(check (M_PI / sqrt(2) * vector3_t (1,1,0)));
  BOOST_CHECK(check (M_PI / sqrt(2) * vector3_t (0,1,1)));
  BOOST_CHECK(check (M_PI / sqrt(2) * vector3_t (1,0,1)));

  BOOST_CHECK(check (M_PI / sqrt(3) * vector3_t (1,-1,1), eps));
  BOOST_CHECK(check (M_PI / sqrt(3) * vector3_t (1, 1,1), eps));

  BOOST_CHECK(check (1e-7 * vector3_t (1,1,1)));
  BOOST_CHECK(check ((M_PI - 0.1 * dlUB) / sqrt(3) * vector3_t (1,-1,1), eps));
  BOOST_CHECK(check ((M_PI - 0.1 * dlUB) / sqrt(3) * vector3_t (1, 1,1), eps));

  BOOST_CHECK(check ((M_PI - 1.001 * dlUB) / sqrt(3) * vector3_t (1, 1,1), eps));
  BOOST_CHECK(check ((M_PI - 0.999 * dlUB) / sqrt(3) * vector3_t (1, 1,1)));
}
