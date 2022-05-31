// Copyright (c) 2020 CNRS, Airbus SAS
// Author: Florent Lamiraux
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr), Florent Lamiraux
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

#ifndef HPP_CONSTRAINTS_COMOARISON_TYPES_HH
#define HPP_CONSTRAINTS_COMOARISON_TYPES_HH

#include <hpp/constraints/fwd.hh>
#include <iostream>
#include <vector>

namespace hpp {
namespace constraints {
namespace internal {
/// \cond HIDDEN_SYMBOLS
struct ReplicateCompType {
  ComparisonType type;
  std::size_t n;

  void push_to(ComparisonTypes_t& v) const {
    for (std::size_t i = 0; i < n; ++i) v.push_back(type);
  }
  // Cast to ComparisonTypes_t
  operator ComparisonTypes_t() const { return ComparisonTypes_t(n, type); }
};  // struct ReplicateCompType
/// \endcond
}  // namespace internal
inline internal::ReplicateCompType operator*(const int& n,
                                             const ComparisonType& c) {
  internal::ReplicateCompType cts;
  cts.type = c;
  cts.n = (std::size_t)n;
  return cts;
}

inline ComparisonTypes_t operator<<(const ComparisonType& a,
                                    const ComparisonType& b) {
  ComparisonTypes_t v(2);
  v[0] = a;
  v[1] = b;
  return v;
}

inline ComparisonTypes_t operator<<(const internal::ReplicateCompType& a,
                                    const ComparisonType& b) {
  ComparisonTypes_t v;
  v.reserve(a.n + 1);
  a.push_to(v);
  v.push_back(b);
  return v;
}

inline ComparisonTypes_t operator<<(const ComparisonType& a,
                                    const internal::ReplicateCompType& b) {
  ComparisonTypes_t v;
  v.reserve(b.n + 1);
  v.push_back(a);
  b.push_to(v);
  return v;
}

inline ComparisonTypes_t& operator<<(ComparisonTypes_t& v,
                                     const ComparisonType& c) {
  v.push_back(c);
  return v;
}

inline ComparisonTypes_t& operator<<(ComparisonTypes_t& v,
                                     const internal::ReplicateCompType& c) {
  for (std::size_t i = 0; i < c.n; ++i) v.push_back(c.type);
  return v;
}

inline ComparisonTypes_t operator<<(const ComparisonTypes_t& v,
                                    const ComparisonType& c) {
  ComparisonTypes_t vv;
  vv.reserve(v.size() + 1);
  vv.insert(vv.end(), v.begin(), v.end());
  vv.push_back(c);
  return vv;
}

inline ComparisonTypes_t operator<<(const ComparisonTypes_t& v,
                                    const internal::ReplicateCompType& c) {
  ComparisonTypes_t vv;
  vv.reserve(v.size() + c.n);
  vv.insert(vv.end(), v.begin(), v.end());
  for (std::size_t i = 0; i < c.n; ++i) vv.push_back(c.type);
  return vv;
}

inline bool operator==(const ComparisonTypes_t& v,
                       const internal::ReplicateCompType& r) {
  if (v.size() != r.n) return false;
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (v[i] != r.type) return false;
  }
  return true;
}

inline std::ostream& operator<<(std::ostream& os,
                                const ComparisonTypes_t& comp) {
  os << "(";
  for (ComparisonTypes_t::const_iterator it = comp.begin(); it != comp.end();
       ++it) {
    if (*it == Equality) {
      os << "Equality";
    } else if (*it == EqualToZero) {
      os << "EqualToZero";
    } else if (*it == Superior) {
      os << "Superior";
    } else if (*it == Inferior) {
      os << "Inferior";
    } else {
      assert("false && ComparisonType out of range.");
    }
    if (it + 1 != comp.end()) os << ", ";
  }
  os << ")";
  return os;
}
}  // namespace constraints
}  // namespace hpp

inline hpp::constraints::ComparisonTypes_t operator<<(
    hpp::constraints::internal::ReplicateCompType c1,
    hpp::constraints::internal::ReplicateCompType c2) {
  hpp::constraints::ComparisonTypes_t vv(c1);
  for (std::size_t i = 0; i < c2.n; ++i) vv.push_back(c2.type);
  return vv;
}

#endif  // HPP_CONSTRAINTS_COMOARISON_TYPES_HH
