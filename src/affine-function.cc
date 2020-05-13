// Copyright (c) 2020, Joseph Mirabel
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

#include <hpp/constraints/affine-function.hh>

#include <pinocchio/serialization/eigen.hpp>

#include <hpp/util/serialization.hh>
#include <hpp/pinocchio/liegroup-element.hh>

namespace hpp {
namespace constraints {
using namespace boost::serialization;

template<class Archive>
void Identity::serialize(Archive & ar, const unsigned int version)
{
  (void) version;
  ar & make_nvp("base", base_object<DifferentiableFunction> (*this));
}

HPP_SERIALIZATION_IMPLEMENT(Identity);

template<class Archive>
void AffineFunction::serialize(Archive & ar, const unsigned int version)
{
  (void) version;
  ar & make_nvp("base", base_object<DifferentiableFunction> (*this));
  ar & make_nvp("J_", const_cast<matrix_t&>(J_));
  ar & make_nvp("b_", const_cast<vector_t&>(b_));
}

HPP_SERIALIZATION_IMPLEMENT(AffineFunction);

template<class Archive>
void ConstantFunction::serialize(Archive & ar, const unsigned int version)
{
  (void) version;
  ar & make_nvp("base", base_object<DifferentiableFunction> (*this));
  ar & make_nvp("c_", const_cast<LiegroupElement&>(c_));
}

HPP_SERIALIZATION_IMPLEMENT(ConstantFunction);

} // namespace constraints
} // namespace hpp

BOOST_CLASS_EXPORT(hpp::constraints::Identity)
BOOST_CLASS_EXPORT(hpp::constraints::AffineFunction)
BOOST_CLASS_EXPORT(hpp::constraints::ConstantFunction)
