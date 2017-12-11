// Copyright (c) 2017 CNRS
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)  and Florent Lamiraux
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

#include <hpp/constraints/matrix-view.hh>

namespace Eigen {
  void BlockIndex::sort (segments_t& a)
  {
    std::sort (a.begin(), a.end(), internal::BlockIndexCompFull ());
  }

  void BlockIndex::shrink (segments_t& a)
  {
    if (a.size() < 2) return;
    // Find consecutive element which overlap
    typename segments_t::iterator e2 = a.begin();
    typename segments_t::iterator e1 = e2++;
    internal::BlockIndexComp<false, true> lend_before_rstart;
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

  bool BlockIndex::overlap (const segment_t& a, const segment_t& b)
  {
    if (a.second == 0) return false;
    if (b.second == 0) return false;

    size_type aend = a.first + a.second;
    size_type bend = b.first + b.second;

    if (a.first <= b.first && b.first < aend) return true;
    if (a.first < bend && bend <= aend)       return true;
    return false;
  }

  size_type BlockIndex::cardinal (const segments_t& a)
  {
    size_type c = 0;
    for (typename segments_t::const_iterator _a = a.begin(); _a != a.end();
         ++_a) c += _a->second;
    return c;
  }

  BlockIndex::segments_t BlockIndex::sum (const segment_t& a,
                                           const segment_t& b)
  {
    if (a.first > b.first) return sum (b, a);
    // a.first <= b.first
    segments_t s (1, a);
    if (a.first + a.second >= b.first)
      s[0].second = std::max (a.second, b.first + b.second - a.first);
    else
      s.push_back (b);
    return s;
  }

  void BlockIndex::add (segments_t& a, const segment_t& b)
  {
    // Sorted insertion of b into a
    segments_t::iterator _it = std::upper_bound (a.begin(), a.end(), b,
        internal::BlockIndexCompFull ());
    a.insert (_it, b);
  }

  void BlockIndex::add (segments_t& a, const segments_t& b)
  {
    // Sorted insertion of b into a, assuming b is sorted.
    segments_t::iterator _a = a.begin();
    for (segments_t::const_iterator _b = b.begin(); _b != b.end(); ++_b) {
      _a = std::upper_bound (_a, a.end(), *_b, internal::BlockIndexCompFull ());
      a.insert (_a, *_b);
    }
  }

  BlockIndex::segments_t BlockIndex::difference (const segment_t& a,
                                                  const segment_t& b)
  {
    if (a.second == 0) return segments_t(0);
    if (b.second == 0) return segments_t(1, a);

    size_type aend = a.first + a.second;
    size_type bend = b.first + b.second;
    segments_t diffs;

    if (a.first < b.first) {
      size_type end = std::min (aend, b.first);
      diffs.push_back (segment_t(a.first, end - a.first));
    }
    if (bend < aend) {
      size_type start = std::max (a.first, bend);
      diffs.push_back (segment_t(start, aend - start));
    }
    return diffs;
  }

  BlockIndex::segments_t BlockIndex::difference (const segments_t& a,
                                                  const segment_t& b)
  {
    segments_t::const_iterator first
      (std::upper_bound (a.begin(), a.end(), b,
                         internal::BlockIndexComp<true, false>()));
    segments_t::const_iterator last
      (std::upper_bound (a.begin(), a.end(), b,
                         internal::BlockIndexComp<false, true>()));
    assert (first == last || last == a.end() || (first != a.end() && first->first + first->second >= last->first));
    segments_t ret; ret.reserve(a.size() + 2);
    ret.insert(ret.end(), a.begin(), first);
    for (typename segments_t::const_iterator _a = first; _a != last; ++_a) {
      segments_t diff = difference (*_a, b);
      ret.insert(ret.end(), diff.begin(), diff.end());
    }
    ret.insert(ret.end(), last, a.end());
    return ret;
  }

  BlockIndex::segments_t BlockIndex::difference (const segment_t& a,
                                                  const segments_t& b)
  {
    segments_t diff (1, a);
    for (typename segments_t::const_iterator _b = b.begin(); _b != b.end();
         ++_b)
      diff = difference (diff, *_b);
    return diff;
  }

  BlockIndex::segments_t BlockIndex::difference (const segments_t& a,
                                                  const segments_t& b)
  {
    segments_t diff;
    for (typename segments_t::const_iterator _a = a.begin(); _a != a.end();
         ++_a) {
      segments_t d = difference(*_a, b);
      diff.insert(diff.end(), d.begin(), d.end());
    }
    return diff;
  }

  BlockIndex::segments_t BlockIndex::split (segments_t& segments,
                                            const size_type& cardinal)
  {
    assert (BlockIndex::cardinal (segments) >= cardinal);
    assert (cardinal >= 0);
    segments_t result;
    size_type remaining = cardinal;
    segments_t::iterator it (segments.begin ());
    while (it != segments.end ()) {
      if (it->second > remaining) {
        result.push_back (segment_t (it->first, remaining));
        it->first += remaining;
        it->second -= remaining;
        return result;
      } else if (it->second == remaining) {
        result.push_back (*it);
        it = segments.erase (it);
        return result;
      } else {
        remaining -= it->second;
        result.push_back (*it);
        it = segments.erase (it);
      }
    }
    return result;
  }

  BlockIndex::segments_t BlockIndex::extract
  (const segments_t& segments, size_type start, size_type cardinal)
  {
    assert (start >= 0);
    assert (cardinal >= 0);
    assert (segments [0].first >= 0);
    assert (start + cardinal <= BlockIndex::cardinal (segments));
    segments_t result;
    size_type remainingToStart = start;
    size_type startIndex;
    size_type remainingToEnd;
    segments_t::const_iterator it (segments.begin ());
    while (it != segments.end ()) {
      if (remainingToStart >= 0) {
        // looking for start
        if (it->second > remainingToStart) {
          startIndex = it->first + remainingToStart;
          if (it->second - remainingToStart >= cardinal) {
            result.push_back (segment_t (startIndex, cardinal));
            return result;
          } else {
            result.push_back (segment_t (startIndex,
                                         it->second - remainingToStart));
            remainingToEnd = cardinal - (it->second - remainingToStart);
            ++it;
            remainingToStart = -1;
          }
        } else {
          remainingToStart -= it->second;
          ++it;
        }
      } else {
        // looking for end
        if (it->second >= remainingToEnd) {
          result.push_back (segment_t (it->first, remainingToEnd));
          return result;
        } else {
          result.push_back (segment_t (it->first, it->second));
          remainingToEnd -= it->second;
          ++it;
        }
      }
    }
    assert (false && "Failed to extract segments_t");
  }
} // namespace Eigen
