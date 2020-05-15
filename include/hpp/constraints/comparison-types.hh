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


#ifndef HPP_CONSTRAINTS_COMPARISON_TYPES_HH
#define HPP_CONSTRAINTS_COMPARISON_TYPES_HH

#include <hpp/constraints/config.hh>
#include <hpp/constraints/fwd.hh>

namespace hpp {
  namespace constraints {
    struct HPP_CONSTRAINTS_DLLAPI ComparisonType
    {
    public:
      ComparisonType (int v) : value(v)
      {
      }
      bool operator==(const ComparisonType& other) const
      {
        return value == other.value;
      }
      bool operator!=(const ComparisonType& other) const
      {
        return value != other.value;
      }
      int value;
    }; // struct ComparisonType
    static ComparisonType Equality(0), EqualToZero(1), Superior(2), Inferior(3);

    /// Vector of ComparisonType.
    ///
    /// This class used to be a vector of ComparisonTypes.
    /// To avoid ambiguity in operator, and operator<<, it is now a class
    /// that contains a vector of comparison types.
    /// Part of the API emulates a vector for compatibility reasons.
    class HPP_CONSTRAINTS_DLLAPI ComparisonTypes
    {
    public:
      typedef std::vector<ComparisonType> Vector_t;
      typedef Vector_t::const_iterator const_iterator;
      typedef Vector_t::iterator iterator;

      /// Return a vector with n times the same value
      /// \param n number of instances,
      /// \param c instance that is copied n times.
      /// \note this method is aimed at replacing a constructor.
      static ComparisonTypes nTimes(std::size_t n, ComparisonType c);
      /// Constructor
      ComparisonTypes();
      /// Constructor with one element
      ComparisonTypes(const ComparisonType& c);
      ComparisonTypes(const ComparisonTypes& other);
      /// \name std::vector API
      /// \{
      std::size_t size() const;
      ComparisonType& operator[](std::size_t i);
      const ComparisonType& operator[](std::size_t i) const;
      bool empty() const;
      void clear();
      void push_back(const ComparisonType& c);
      iterator insert(const_iterator pos, const_iterator first,
                      const_iterator last );
      iterator begin();
      const_iterator begin() const;
      iterator end();
      const_iterator end() const;
      friend bool operator==(const ComparisonTypes& lhs,
                             const ComparisonTypes& rhs);
      friend bool operator!=(const ComparisonTypes& lhs,
                             const ComparisonTypes& rhs);
      /// \}
    private:
      Vector_t vector_;
    }; // class ComparisonTypes
    bool operator==(const ComparisonTypes& lhs, const ComparisonTypes& rhs);
    bool operator!=(const ComparisonTypes& lhs, const ComparisonTypes& rhs);
    ComparisonTypes operator<<(ComparisonType c1, ComparisonType c2);
    ComparisonTypes operator<<(const ComparisonTypes& comp1,
                               const ComparisonTypes& comp2);
    ComparisonTypes operator<<(const ComparisonTypes& comp, ComparisonType c);
    ComparisonTypes operator*(std::size_t n, ComparisonType c);
    std::ostream& operator<<(std::ostream& os, const ComparisonTypes& comp);
  } // namespace constraints
} // namespace hpp
#endif //HPP_CONSTRAINTS_COMPARISON_TYPES_HH
