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

#ifndef HPP_CONSTRAINTS_SERIALIZATION_HH
#define HPP_CONSTRAINTS_SERIALIZATION_HH

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

#include <hpp/constraints/matrix-view.hh>

namespace boost {
namespace serialization {

template<class Archive> inline void segments_serialize(Archive & ar, const char* n, Eigen::BlockIndex::segments_t& s) { ar & make_nvp(n, s); }
template<class Archive> inline void segments_serialize(Archive &, const char*, Eigen::internal::empty_struct&) {}

template<class Archive, bool _allRows, bool _allCols>
inline void serialize(Archive & ar, Eigen::MatrixBlocks<_allRows, _allCols>& b,
    const unsigned int version)
{
  (void) version;
  ar & make_nvp("nbRows", b.m_nbRows);
  ar & make_nvp("nbCols", b.m_nbCols);
  segments_serialize(ar, "rows", b.m_rows);
  segments_serialize(ar, "cols", b.m_cols);
}
} // namespace serialization
} // namespace boost

#endif // HPP_CONSTRAINTS_SERIALIZATION_HH
