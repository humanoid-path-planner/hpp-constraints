// Copyright (c) 2020, CNRS - Airbus SAS
// Authors: Florent Lamiraux
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

#include <hpp/pinocchio/liegroup-space.hh>
#include <hpp/constraints/affine-function.hh>
#include <hpp/constraints/implicit.hh>

#define BOOST_TEST_MODULE ComparisonTypes
#include <boost/test/included/unit_test.hpp>

using hpp::constraints::Implicit;
using hpp::constraints::ImplicitPtr_t;
using hpp::constraints::DifferentiableFunctionPtr_t;
using hpp::constraints::ComparisonType;
using hpp::constraints::ComparisonTypes_t;
using hpp::constraints::Equality;
using hpp::constraints::EqualToZero;
using hpp::constraints::Superior;
using hpp::constraints::Inferior;
using hpp::constraints::Identity;
using hpp::pinocchio::LiegroupSpace;

BOOST_AUTO_TEST_CASE (operators)
//int main()
{
  // [Inferior]
  ComparisonTypes_t expected, actual(1, Inferior);
  expected.push_back(Inferior);
  BOOST_CHECK_EQUAL(expected, actual);
  // [Equality, EqualToZero]
  expected.clear(); actual.clear();
  expected.push_back(Equality);
  expected.push_back(EqualToZero);
  actual = ComparisonTypes_t(Equality << EqualToZero);
  BOOST_CHECK_EQUAL(expected, actual);
  BOOST_CHECK_EQUAL(expected, Equality << EqualToZero);

  // [Equality, EqualToZero, Inferior]
  expected.clear(); actual.clear();
  expected.push_back(Equality);
  expected.push_back(EqualToZero);
  expected.push_back(Inferior);
  actual = (Equality << EqualToZero << Inferior);
  BOOST_CHECK_EQUAL(expected, actual);
  BOOST_CHECK_EQUAL(expected, Equality << EqualToZero << Inferior);

  // [Equality, EqualToZero, 4 * Inferior]
  expected.clear(); actual.clear();
  expected.push_back(Equality);
  expected.push_back(EqualToZero);
  expected.push_back(Inferior);
  expected.push_back(Inferior);
  expected.push_back(Inferior);
  expected.push_back(Inferior);

  actual = (Equality << EqualToZero << 4 * Inferior);
  BOOST_CHECK_EQUAL(expected, actual);
  BOOST_CHECK_EQUAL(expected, Equality << EqualToZero << 4 * Inferior);
  // [EqualToZero, EqualToZero, EqualToZero, EqualToZero, EqualToZero,
  //  Equality, Equality, Equality] -> ConvexShapeContact
  expected.clear(); actual.clear();
  expected.push_back(EqualToZero);
  expected.push_back(EqualToZero);
  expected.push_back(EqualToZero);
  expected.push_back(EqualToZero);
  expected.push_back(EqualToZero);
  expected.push_back(Equality);
  expected.push_back(Equality);
  expected.push_back(Equality);

  actual = (EqualToZero << EqualToZero << EqualToZero << EqualToZero
            << EqualToZero << Equality << Equality << Equality);
  BOOST_CHECK_EQUAL(expected, actual);
  ImplicitPtr_t constraint(Implicit::create(Identity::create
                                            (LiegroupSpace::Rn(8), "I8"),
    EqualToZero << EqualToZero << EqualToZero << EqualToZero << EqualToZero
    << Equality << Equality << Equality));
  BOOST_CHECK_EQUAL(expected, constraint->comparisonType());
  actual = ComparisonTypes_t(EqualToZero << EqualToZero << EqualToZero
                             << EqualToZero << EqualToZero << Equality
                             << Equality << Equality);
  BOOST_CHECK_EQUAL(expected, actual);
}
