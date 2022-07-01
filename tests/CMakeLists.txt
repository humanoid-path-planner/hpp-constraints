# Copyright 2012, 2013, 2014 CNRS-LAAS
#
# Author: Joseph Mirabel
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

# Make Boost.Test generates the main function in test cases.
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

# ADD_TESTCASE(NAME)
# ------------------------
#
# Define a test named `NAME'.
#
# This macro will create a binary from `NAME.cc', link it against
# Boost and add it to the test suite.
#

find_library(RT_LIB rt)

MACRO(ADD_TESTCASE NAME)
  ADD_UNIT_TEST(${NAME} ${NAME}.cc)
  TARGET_LINK_LIBRARIES(${NAME} PRIVATE ${PROJECT_NAME} Boost::unit_test_framework)
  if(RT_LIB)
    TARGET_LINK_LIBRARIES(${NAME} PRIVATE ${RT_LIB})
  endif()
ENDMACRO(ADD_TESTCASE)

FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ADD_TESTCASE(multithread)
ENDIF()

ADD_TESTCASE(comparison-types)
ADD_TESTCASE(logarithm)
ADD_TESTCASE(matrix-view)
ADD_TESTCASE(svd)

ADD_TESTCASE(test-jacobians)
ADD_TESTCASE(generic-transformation)
ADD_TESTCASE(convex-shape)
ADD_TESTCASE(convex-shape-contact)
ADD_TESTCASE(symbolic-calculus)
ADD_TESTCASE(solver-hierarchical-iterative)
ADD_TESTCASE(explicit-constraint-set)
ADD_TESTCASE(solver-by-substitution)
ADD_TESTCASE(gjk)
