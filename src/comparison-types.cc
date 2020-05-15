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

#include <hpp/constraints/comparison-types.hh>

namespace hpp {
  namespace constraints {
    //ComparisonType Equality(0), EqualToZero(1), Superior(2), Inferior(3);

    ComparisonTypes ComparisonTypes::nTimes(std::size_t n, ComparisonType c)
    {
      ComparisonTypes res;
      res.vector_ = Vector_t(n, c);
      return res;
    }

    ComparisonTypes::ComparisonTypes(): vector_()
    {
    }
    ComparisonTypes::ComparisonTypes(const ComparisonType& c): vector_()
    {
      vector_.push_back(c);
    }
    ComparisonTypes::ComparisonTypes(const ComparisonTypes& other):
      vector_(other.vector_)
    {
    }
    std::size_t ComparisonTypes::size() const
    {
      return vector_.size();
    }
    ComparisonType& ComparisonTypes::operator[](std::size_t i)
    {
      return vector_[i];
    }
    const ComparisonType& ComparisonTypes::operator[](std::size_t i) const
    {
      return vector_[i];
    }
    bool ComparisonTypes::empty() const
    {
      return vector_.empty();
    }
    void ComparisonTypes::clear()
    {
      vector_.clear();
    }
    void ComparisonTypes::push_back(const ComparisonType& c)
    {
      vector_.push_back(c);
    }

    ComparisonTypes::iterator ComparisonTypes::insert
    (const_iterator pos, const_iterator first, const_iterator last )
    {
      return vector_.insert(pos, first, last);
    }
    
    ComparisonTypes::iterator ComparisonTypes::begin()
    {
      return vector_.begin();
    }
    ComparisonTypes::const_iterator ComparisonTypes::begin() const
    {
      return vector_.begin();
    }
    ComparisonTypes::iterator ComparisonTypes::end()
    {
      return vector_.end();
    }
    ComparisonTypes::const_iterator ComparisonTypes::end() const
    {
      return vector_.end();
    }
 
    bool operator!=(const ComparisonTypes& lhs,
                    const ComparisonTypes& rhs)
    {
      return (lhs.vector_ != rhs.vector_);
    }

    bool operator==(const ComparisonTypes& lhs,
                    const ComparisonTypes& rhs)
    {
      return (lhs.vector_ == rhs.vector_);
    }
    // Assignment of comparison types
    ComparisonTypes operator<<(ComparisonType c1, ComparisonType c2)
    {
      ComparisonTypes res;
      res.push_back(c1);
      res.push_back(c2);
      return res;
    }

    ComparisonTypes operator<<(const ComparisonTypes& comp1,
                               const ComparisonTypes& comp2)
    {
      ComparisonTypes res(comp1);
      res.insert(res.end(), comp2.begin(), comp2.end());
      return res;
    }

    ComparisonTypes operator<<(const ComparisonTypes& comp, ComparisonType c)
    {
      ComparisonTypes res(comp);
      res.push_back(c);
      return res;
    }

    ComparisonTypes operator*(std::size_t n, ComparisonType c)
    {
      return ComparisonTypes::nTimes(n, c);
    }
    std::ostream& operator<<(std::ostream& os, const ComparisonTypes& comp)
    {
      os << "(";
      for (ComparisonTypes::const_iterator it=comp.begin();
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
