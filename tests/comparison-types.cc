// Copyright (c) 2020, CNRS - Airbus SAS
// Authors: Florent Lamiraux
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

#include <hpp/constraints/affine-function.hh>
#include <hpp/constraints/implicit.hh>
#include <hpp/pinocchio/liegroup-space.hh>

#define BOOST_TEST_MODULE ComparisonTypes
#include <boost/test/included/unit_test.hpp>

using hpp::constraints::ComparisonType;
using hpp::constraints::ComparisonTypes_t;
using hpp::constraints::DifferentiableFunctionPtr_t;
using hpp::constraints::Equality;
using hpp::constraints::EqualToZero;
using hpp::constraints::Identity;
using hpp::constraints::Implicit;
using hpp::constraints::ImplicitPtr_t;
using hpp::constraints::Inferior;
using hpp::constraints::Superior;
using hpp::pinocchio::LiegroupSpace;

BOOST_AUTO_TEST_CASE(operators)
// int main()
{
  // [Inferior]
  ComparisonTypes_t expected, actual(1, Inferior);
  expected.push_back(Inferior);
  BOOST_CHECK_EQUAL(expected, actual);
  // [Equality, EqualToZero]
  expected.clear();
  actual.clear();
  expected.push_back(Equality);
  expected.push_back(EqualToZero);
  actual = ComparisonTypes_t(Equality << EqualToZero);
  BOOST_CHECK_EQUAL(expected, actual);
  BOOST_CHECK_EQUAL(expected, Equality << EqualToZero);

  // [Equality, EqualToZero, Inferior]
  expected.clear();
  actual.clear();
  expected.push_back(Equality);
  expected.push_back(EqualToZero);
  expected.push_back(Inferior);
  actual = (Equality << EqualToZero << Inferior);
  BOOST_CHECK_EQUAL(expected, actual);
  BOOST_CHECK_EQUAL(expected, Equality << EqualToZero << Inferior);

  // [Equality, EqualToZero, 4 * Inferior]
  expected.clear();
  actual.clear();
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
  expected.clear();
  actual.clear();
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
  ImplicitPtr_t constraint(Implicit::create(
      Identity::create(LiegroupSpace::Rn(8), "I8"),
      EqualToZero << EqualToZero << EqualToZero << EqualToZero << EqualToZero
                  << Equality << Equality << Equality));
  BOOST_CHECK_EQUAL(expected, constraint->comparisonType());
  actual = ComparisonTypes_t(EqualToZero << EqualToZero << EqualToZero
                                         << EqualToZero << EqualToZero
                                         << Equality << Equality << Equality);
  BOOST_CHECK_EQUAL(expected, actual);
}
