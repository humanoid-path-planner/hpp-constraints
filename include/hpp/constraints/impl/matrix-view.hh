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
  template <typename Derived>
  typename BlockIndex::segments_t BlockIndex::fromLogicalExpression
  (const Eigen::ArrayBase<Derived>& array)
  {
    segments_t res;
    for (size_type i = 0; i < array.derived().size(); ++i)
      if (array.derived()[i]) res.push_back (segment_t(i, 1));
    shrink(res);
    return res;
  }
} // namespace Eigen

#include <hpp/pinocchio/util.hh>

namespace hpp {
  template <int Option> struct prettyPrint <Eigen::BlockIndex::segments_t, Option>
  {
    static std::ostream& run (std::ostream& os, const constraints::segments_t& segs)
    {
      Eigen::internal::print_indices::run (os, segs);
      return os;
    }
  };
} // namespace hpp
