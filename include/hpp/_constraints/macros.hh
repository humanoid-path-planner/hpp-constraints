// Copyright (c) 2015 CNRS
// Author: Joseph Mirabel
//
//
// This file is part of hpp-constraints
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
// hpp-constraints  If not, see
// <http://www.gnu.org/licenses/>.

#ifndef HPP_CONSTRAINTS_MACROS_HH
# define HPP_CONSTRAINTS_MACROS_HH

# ifdef HPP_DEBUG

#  define HPP_DEBUG_SVDCHECK(svd)                                       \
  do {                                                                  \
    if (svd.rank () > 0) {                                              \
      value_type SSV = svd.singularValues()(svd.rank()-1);              \
      if (std::abs (SSV) < 1e-8) {                                      \
        hppDout (warning, "SVD check - low singular value: " << SSV);   \
      }                                                                 \
    }                                                                   \
  } while (0)

#  ifdef HPP_CONSTRAINTS_NUMERIC_DEBUG

#   define hppDnum(channel, data) hppDout (channel, data)

#  else // HPP_CONSTRAINTS_NUMERIC_DEBUG

#   define hppDnum(channel, data) do { } while (0)

#  endif // HPP_CONSTRAINTS_NUMERIC_DEBUG

# else // HPP_DEBUG

#  define HPP_DEBUG_SVDCHECK(svd) do { } while (0)
#  define hppDnum(channel, data) do { } while (0)

# endif // HPP_DEBUG

#endif // HPP_CONSTRAINTS_MACROS_HH
