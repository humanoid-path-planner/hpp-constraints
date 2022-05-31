// Copyright (c) 2020, Joseph Mirabel
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

#ifndef HPP_CONSTRAINTS_SERIALIZATION_HH
#define HPP_CONSTRAINTS_SERIALIZATION_HH

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <hpp/constraints/matrix-view.hh>

namespace boost {
namespace serialization {

template <class Archive>
inline void segments_serialize(Archive& ar, const char* n,
                               Eigen::BlockIndex::segments_t& s) {
  ar& make_nvp(n, s);
}
template <class Archive>
inline void segments_serialize(Archive&, const char*,
                               Eigen::internal::empty_struct&) {}

template <class Archive, bool _allRows, bool _allCols>
inline void serialize(Archive& ar, Eigen::MatrixBlocks<_allRows, _allCols>& b,
                      const unsigned int version) {
  (void)version;
  ar& make_nvp("nbRows", b.m_nbRows);
  ar& make_nvp("nbCols", b.m_nbCols);
  segments_serialize(ar, "rows", b.m_rows);
  segments_serialize(ar, "cols", b.m_cols);
}
}  // namespace serialization
}  // namespace boost

#endif  // HPP_CONSTRAINTS_SERIALIZATION_HH
