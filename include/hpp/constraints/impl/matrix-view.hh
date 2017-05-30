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
  template <typename IndexType>
  bool BlockIndex<IndexType>::overlap (const type& a, const type& b)
  {
    if (a.second == 0) return false;
    if (b.second == 0) return false;

    Index aend = a.first + a.second;
    Index bend = b.first + b.second;

    if (a.first <= b.first && b.first < aend) return true;
    if (a.first < bend && bend <= aend)       return true;
    return false;
  }

  template <typename IndexType>
  typename BlockIndex<IndexType>::vector_t BlockIndex<IndexType>::difference (const type& a, const type& b)
  {
    if (a.second == 0) return vector_t(0);
    if (b.second == 0) return vector_t(1, a);

    Index aend = a.first + a.second;
    Index bend = b.first + b.second;
    vector_t diffs;

    if (a.first < b.first) {
      Index end = std::min (aend, b.first);
      diffs.push_back (type(a.first, end - a.first));
    }
    if (bend < aend) {
      Index start = std::max (a.first, bend);
      diffs.push_back (type(start, aend - start));
    }
    return diffs;
  }
} // namespace Eigen
