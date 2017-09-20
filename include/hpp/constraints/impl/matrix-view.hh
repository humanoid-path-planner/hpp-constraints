// Copyright (c) 2017, Joseph Mirabel
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

namespace Eigen {

  typedef hpp::constraints::size_type size_type;

  namespace internal {
    template <bool lfirst, bool rfirst>
    struct BlockIndexComp {
      typedef BlockIndex::interval_t interval_t;
      bool operator() (const interval_t& l, const interval_t& r) const {
        return ( lfirst ? l.first : l.first + l.second )
          <    ( rfirst ? r.first : r.first + r.second );
      }
    };
    struct BlockIndexCompFull {
      typedef BlockIndex::interval_t interval_t;
      bool operator() (const interval_t& l, const interval_t& r) const {
        return ( l.first  < r.first )
          ||   ( l.first == r.first && l.second < r.second );
      }
    };
  }

  template <typename Derived>
  typename BlockIndex::vector_t BlockIndex::fromLogicalExpression
  (const Eigen::ArrayBase<Derived>& array)
  {
    vector_t res;
    for (std::size_t i = 0; i < array.derived().size(); ++i)
      if (array.derived()[i]) res.push_back (interval_t(i, 1));
    shrink(res);
    return res;
  }

} // namespace Eigen
