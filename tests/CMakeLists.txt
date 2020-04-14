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
MACRO(ADD_TESTCASE NAME)
  ADD_UNIT_TEST(${NAME} ${NAME}.cc)
  TARGET_LINK_LIBRARIES(${NAME} ${Boost_LIBRARIES} ${PROJECT_NAME})
ENDMACRO(ADD_TESTCASE)

FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ADD_TESTCASE(multithread)
ENDIF()

ADD_TESTCASE(logarithm)
ADD_TESTCASE(matrix-view)
ADD_TESTCASE(svd)

ADD_TESTCASE(test-jacobians)
ADD_TESTCASE(generic-transformation)
ADD_TESTCASE(convex-shape)
ADD_TESTCASE(symbolic-calculus)
ADD_TESTCASE(solver-hierarchical-iterative)
ADD_TESTCASE(explicit-constraint-set)
ADD_TESTCASE(solver-by-substitution)
ADD_TESTCASE(gjk)
