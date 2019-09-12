# Copyright 2012, 2013, 2014 CNRS-LAAS
#
# Author: Joseph Mirabel
#
# This file is part of hpp-constraints
# hpp-constraints is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hpp-constraints is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Lesser Public License for more details.
# You should have received a copy of the GNU Lesser General Public License
# along with hpp-constraints  If not, see <http://www.gnu.org/licenses/>.

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

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
MACRO(ADD_TESTCASE NAME GENERATED HEADER_ONLY)
  IF (${GENERATED} STREQUAL TRUE)
    ADD_UNIT_TEST(${NAME} ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.cc)
  ELSE()
    ADD_UNIT_TEST(${NAME} ${NAME}.cc)
  ENDIF()

  PKG_CONFIG_USE_DEPENDENCY(${NAME} hpp-pinocchio)

  # Link against Boost and project library.
  IF(${HEADER_ONLY})
    TARGET_LINK_LIBRARIES(${NAME} ${Boost_LIBRARIES})
  ELSE(${HEADER_ONLY})
    TARGET_LINK_LIBRARIES(${NAME} ${Boost_LIBRARIES} hpp-constraints)
  ENDIF(${HEADER_ONLY})
ENDMACRO(ADD_TESTCASE)

FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ADD_TESTCASE (multithread             FALSE FALSE)
ENDIF()

ADD_TESTCASE (logarithm                 FALSE TRUE )
ADD_TESTCASE (matrix-view               FALSE FALSE )
ADD_TESTCASE (svd                       FALSE TRUE )

ADD_TESTCASE (test-jacobians            FALSE FALSE)
ADD_TESTCASE (generic-transformation    FALSE FALSE)
ADD_TESTCASE (convex-shape              FALSE FALSE)
ADD_TESTCASE (symbolic-calculus         FALSE FALSE)
ADD_TESTCASE (solver-hierarchical-iterative FALSE FALSE)
ADD_TESTCASE (explicit-constraint-set   FALSE FALSE)
ADD_TESTCASE (solver-by-substitution    FALSE FALSE)
ADD_TESTCASE (gjk                       FALSE FALSE)
