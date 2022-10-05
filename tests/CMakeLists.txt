# Copyright 2012, 2013, 2014 CNRS-LAAS
#
# Author: Joseph Mirabel
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# 1. Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Make Boost.Test generates the main function in test cases.
add_definitions(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

# ADD_TESTCASE(NAME)
# ------------------------
#
# Define a test named `NAME'.
#
# This macro will create a binary from `NAME.cc', link it against Boost and add
# it to the test suite.
#

find_library(RT_LIB rt)

macro(ADD_TESTCASE NAME)
  add_unit_test(${NAME} ${NAME}.cc)
  target_link_libraries(${NAME} PRIVATE ${PROJECT_NAME}
                                        Boost::unit_test_framework)
  if(RT_LIB)
    target_link_libraries(${NAME} PRIVATE ${RT_LIB})
  endif()
endmacro(ADD_TESTCASE)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  add_testcase(multithread)
endif()

add_testcase(comparison-types)
add_testcase(logarithm)
add_testcase(matrix-view)
add_testcase(svd)

add_testcase(test-jacobians)
add_testcase(generic-transformation)
add_testcase(convex-shape)
add_testcase(convex-shape-contact)
add_testcase(symbolic-calculus)
add_testcase(solver-hierarchical-iterative)
add_testcase(explicit-constraint-set)
add_testcase(solver-by-substitution)
add_testcase(gjk)
