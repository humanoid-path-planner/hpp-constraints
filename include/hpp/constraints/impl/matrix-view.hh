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
  namespace internal {
    template <typename IndexType, bool lfirst, bool rfirst>
    struct BlockIndexComp {
      typedef typename BlockIndex<IndexType>::type type;
      bool operator() (const type& l, const type& r) const {
        return ( lfirst ? l.first : l.first + l.second )
          <    ( rfirst ? r.first : r.first + r.second );
      }
    };
    template <typename IndexType>
    struct BlockIndexCompFull {
      typedef typename BlockIndex<IndexType>::type type;
      bool operator() (const type& l, const type& r) const {
        return ( l.first  < r.first )
          ||   ( l.first == r.first && l.second < r.second );
      }
    };
  }

  template <typename IndexType>
  template <typename Derived>
  typename BlockIndex<IndexType>::vector_t BlockIndex<IndexType>::fromLogicalExpression
  (const Eigen::ArrayBase<Derived>& array)
  {
    vector_t res;
    for (std::size_t i = 0; i < array.derived().size(); ++i)
      if (array.derived()[i]) res.push_back (type(i, 1));
    shrink(res);
    return res;
  }

  template <typename IndexType>
  void BlockIndex<IndexType>::sort (vector_t& a)
  {
    std::sort (a.begin(), a.end(), internal::BlockIndexCompFull<IndexType>());
  }

  template <typename IndexType>
  void BlockIndex<IndexType>::shrink (vector_t& a)
  {
    if (a.size() < 2) return;
    // Find consecutive element which overlap
    typename vector_t::iterator e2 = a.begin();
    typename vector_t::iterator e1 = e2++;
    internal::BlockIndexComp<IndexType, false, true> lend_before_rstart;
    while (e2 != a.end()) {
      if (!lend_before_rstart(*e1, *e2)) {
        e1->second = std::max(e1->second, e2->first + e2->second - e1->first);
        e2 = a.erase(e2);
      } else {
        e1 = e2;
        ++e2;
      }
    }
  }

  template <typename IndexType>
  inline bool BlockIndex<IndexType>::overlap (const type& a, const type& b)
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
  inline IndexType BlockIndex<IndexType>::cardinal (const vector_t& a)
  {
    IndexType c = 0;
    for (typename vector_t::const_iterator _a = a.begin(); _a != a.end(); ++_a) c += _a->second;
    return c;
  }

  template <typename IndexType>
  inline typename BlockIndex<IndexType>::vector_t BlockIndex<IndexType>::sum (const type& a, const type& b)
  {
    if (a.first > b.first) return sum (b, a);
    // a.first <= b.first
    vector_t s (1, a);
    if (a.first + a.second >= b.first)
      s[0].second = std::max (a.second, b.first + b.second - a.first);
    else
      s.push_back (b);
    return s;
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

  template <typename IndexType>
  typename BlockIndex<IndexType>::vector_t BlockIndex<IndexType>::difference (const vector_t& a, const type& b)
  {
    typename vector_t::const_iterator first = std::upper_bound (a.begin(), a.end(), b, internal::BlockIndexComp<IndexType, true, false>());
    typename vector_t::const_iterator last  = std::upper_bound (a.begin(), a.end(), b, internal::BlockIndexComp<IndexType, false, true>());
    assert (first == last || last == a.end() || (first != a.end() && first->first + first->second >= last->first));
    vector_t ret; ret.reserve(a.size() + 2);
    ret.insert(ret.end(), a.begin(), first);
    for (typename vector_t::const_iterator _a = first; _a != last; ++_a) {
      vector_t diff = difference (*_a, b);
      ret.insert(ret.end(), diff.begin(), diff.end());
    }
    ret.insert(ret.end(), last, a.end());
    return ret;
  }

  template <typename IndexType>
  typename BlockIndex<IndexType>::vector_t BlockIndex<IndexType>::difference (const type& a, const vector_t& b)
  {
    vector_t diff (1, a);
    for (typename vector_t::const_iterator _b = b.begin(); _b != b.end(); ++_b)
      diff = difference (diff, *_b);
    return diff;
  }

  template <typename IndexType>
  typename BlockIndex<IndexType>::vector_t BlockIndex<IndexType>::difference (const vector_t& a, const vector_t& b)
  {
    vector_t diff;
    for (typename vector_t::const_iterator _a = a.begin(); _a != a.end(); ++_a) {
      vector_t d = difference(*_a, b);
      diff.insert(diff.end(), d.begin(), d.end());
    }
    return diff;
  }
} // namespace Eigen
