// Copyright (c) 2020 CNRS, Airbus SAS
// Author: Florent Lamiraux
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr), Florent Lamiraux
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

#ifndef HPP_CONSTRAINTS_COMOARISON_TYPES_HH
#define HPP_CONSTRAINTS_COMOARISON_TYPES_HH

#include <iostream>
#include <vector>
#include <hpp/constraints/fwd.hh>

namespace hpp {
  namespace constraints {
    namespace internal {
      /// \cond HIDDEN_SYMBOLS
      struct ReplicateCompType {
        ComparisonType type;
        std::size_t n;

        void push_to(ComparisonTypes_t& v) const
        {
          for (std::size_t i = 0; i < n; ++i) v.push_back(type);
        }
        // Cast to ComparisonTypes_t
        operator ComparisonTypes_t() const
        {
          return ComparisonTypes_t(n, type);
        }
      }; // struct ReplicateCompType
      /// \endcond
    }
    inline internal::ReplicateCompType operator*(const int& n,
                                                 const ComparisonType& c)
    {
      internal::ReplicateCompType cts;
      cts.type = c;
      cts.n = (std::size_t)n;
      return cts;
    }

    inline ComparisonTypes_t operator<<(const ComparisonType& a,
                                        const ComparisonType& b)
    {
      ComparisonTypes_t v(2);
      v[0]=a;
      v[1]=b;
      return v;
    }

    inline ComparisonTypes_t operator<<(const internal::ReplicateCompType& a,
                                        const ComparisonType& b)
    {
      ComparisonTypes_t v;
      v.reserve(a.n+1);
      a.push_to(v);
      v.push_back(b);
      return v;
    }

    inline ComparisonTypes_t operator<<(const ComparisonType& a,
                                        const internal::ReplicateCompType& b)
    {
      ComparisonTypes_t v;
      v.reserve(b.n+1);
      v.push_back(a);
      b.push_to(v);
      return v;
    }

    inline ComparisonTypes_t& operator<<(ComparisonTypes_t& v,
                                         const ComparisonType& c)
    {
      v.push_back(c);
      return v;
    }

    inline ComparisonTypes_t& operator<<(ComparisonTypes_t& v,
                                         const internal::ReplicateCompType& c)
    {
      for (std::size_t i = 0; i < c.n; ++i) v.push_back(c.type);
      return v;
    }

    inline ComparisonTypes_t operator<<(const ComparisonTypes_t& v,
                                        const ComparisonType& c)
    {
      ComparisonTypes_t vv;
      vv.reserve(v.size()+1);
      vv.insert(vv.end(), v.begin(), v.end());
      vv.push_back(c);
      return vv;
    }

    inline ComparisonTypes_t operator<<(const ComparisonTypes_t& v,
                                        const internal::ReplicateCompType& c)
    {
      ComparisonTypes_t vv;
      vv.reserve(v.size()+c.n);
      vv.insert(vv.end(), v.begin(), v.end());
      for (std::size_t i = 0; i < c.n; ++i) vv.push_back(c.type);
      return vv;
    }

    inline bool operator==(const ComparisonTypes_t& v,
                           const internal::ReplicateCompType& r)
    {
      if (v.size() != r.n) return false;
      for (std::size_t i=0; i<v.size(); ++i)
      {
        if (v[i] != r.type) return false;
      }
      return true;
    }

    inline std::ostream& operator<<(std::ostream& os,
                                    const ComparisonTypes_t& comp)
    {
      os << "(";
      for (ComparisonTypes_t::const_iterator it=comp.begin();
           it!=comp.end(); ++it)
      {
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
        if(it+1 != comp.end())
          os << ", ";
      }
      os << ")";
      return os;
    }
  } // namespace constraints
} // namespace hpp

inline hpp::constraints::ComparisonTypes_t operator<<
(hpp::constraints::internal::ReplicateCompType c1,
 hpp::constraints::internal::ReplicateCompType c2)
{
  hpp::constraints::ComparisonTypes_t vv(c1);
  for (std::size_t i = 0; i < c2.n; ++i) vv.push_back(c2.type);
  return vv;
}


#endif // HPP_CONSTRAINTS_COMOARISON_TYPES_HH
