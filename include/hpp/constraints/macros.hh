// Copyright (c) 2015 CNRS
// Author: Joseph Mirabel
//
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

#ifndef HPP_CONSTRAINTS_MACROS_HH
#define HPP_CONSTRAINTS_MACROS_HH

#ifdef HPP_DEBUG

#define HPP_DEBUG_SVDCHECK(svd)                                      \
  do {                                                               \
    if (svd.rank() > 0) {                                            \
      value_type SSV = svd.singularValues()(svd.rank() - 1);         \
      if (std::abs(SSV) < 1e-8) {                                    \
        hppDout(warning, "SVD check - low singular value: " << SSV); \
      }                                                              \
    }                                                                \
  } while (0)

#ifdef HPP_CONSTRAINTS_NUMERIC_DEBUG

#define hppDnum(channel, data) hppDout(channel, data)

#else  // HPP_CONSTRAINTS_NUMERIC_DEBUG

#define hppDnum(channel, data) \
  do {                         \
  } while (0)

#endif  // HPP_CONSTRAINTS_NUMERIC_DEBUG

#else  // HPP_DEBUG

#define HPP_DEBUG_SVDCHECK(svd) \
  do {                          \
  } while (0)
#define hppDnum(channel, data) \
  do {                         \
  } while (0)

#endif  // HPP_DEBUG

#endif  // HPP_CONSTRAINTS_MACROS_HH
