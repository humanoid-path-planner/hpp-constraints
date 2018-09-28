# Copyright 2014, Florent Lamiraux, CNRS
#
# This file is part of hpp-constraints.
# hpp-constraints is free software: you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# hpp-constraints is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Lesser Public License for more details. You should have
# received a copy of the GNU Lesser General Public License along with
# hpp-constraints. If not, see <http://www.gnu.org/licenses/>.

# Requires at least CMake 2.6 to configure the package.
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

CMAKE_POLICY (SET CMP0012 NEW)
SET(CXX_DISABLE_WERROR TRUE)
SET(DOXYGEN_USE_MATHJAX YES)
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/eigen.cmake)
INCLUDE(cmake/test.cmake)

SET(PROJECT_NAME hpp-constraints)
SET(PROJECT_DESCRIPTION
  "Definition of basic geometric constraints for motion planning")
SET(PROJECT_URL "https://github.com/humanoid-path-planner/hpp-constraints")

# Add a cache variable to allow not compiling and running tests
set (RUN_TESTS TRUE CACHE BOOL "compile and run unit tests")

SET(${PROJECT_NAME}_HEADERS
  include/hpp/constraints/differentiable-function.hh
  include/hpp/constraints/differentiable-function-set.hh
  include/hpp/constraints/differentiable-function-stack.hh
  include/hpp/constraints/active-set-differentiable-function.hh
  include/hpp/constraints/affine-function.hh
  include/hpp/constraints/distance-between-bodies.hh
  include/hpp/constraints/fwd.hh
  include/hpp/constraints/svd.hh
  include/hpp/constraints/tools.hh
  include/hpp/constraints/matrix-view.hh
  include/hpp/constraints/generic-transformation.hh
  include/hpp/constraints/macros.hh
  include/hpp/constraints/convex-shape.hh
  include/hpp/constraints/convex-shape-contact.hh
  include/hpp/constraints/symbolic-calculus.hh
  include/hpp/constraints/symbolic-function.hh
  include/hpp/constraints/relative-com.hh
  include/hpp/constraints/com-between-feet.hh
  include/hpp/constraints/configuration-constraint.hh
  include/hpp/constraints/explicit-constraint-set.hh
  include/hpp/constraints/explicit-solver.hh
  include/hpp/constraints/explicit.hh
  include/hpp/constraints/explicit/implicit-function.hh
  include/hpp/constraints/explicit/relative-pose.hh
  include/hpp/constraints/explicit/relative-transformation.hh
  include/hpp/constraints/implicit.hh
  include/hpp/constraints/implicit/relative-pose.hh
  include/hpp/constraints/implicit-constraint-set.hh
  include/hpp/constraints/hybrid-solver.hh
  include/hpp/constraints/iterative-solver.hh
  include/hpp/constraints/locked-joint.hh
  include/hpp/constraints/solver/hierarchical-iterative.hh
  include/hpp/constraints/solver/by-substitution.hh

  include/hpp/constraints/function/of-parameter-subset.hh

  include/hpp/constraints/solver/impl/by-substitution.hh
  include/hpp/constraints/solver/impl/hierarchical-iterative.hh
  include/hpp/constraints/impl/matrix-view.hh
  include/hpp/constraints/impl/matrix-view-operation.hh
)

# Activate hpp-util logging if requested
SET (HPP_DEBUG FALSE CACHE BOOL "trigger hpp-util debug output")
IF (HPP_DEBUG)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHPP_DEBUG")
ENDIF()

# Activate hpp-util logging if requested
SET (HPP_CONSTRAINTS_NUMERIC_DEBUG ${HPP_DEBUG} CACHE BOOL "trigger hpp-util debug output")
IF (HPP_CONSTRAINTS_NUMERIC_DEBUG)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHPP_CONSTRAINTS_NUMERIC_DEBUG")
ENDIF()

OPTION(CHECK_JACOBIANS "Check jacobians against finite difference." OFF)
IF(CHECK_JACOBIANS)
  ADD_DEFINITIONS(-DCHECK_JACOBIANS)
ENDIF(CHECK_JACOBIANS)

# Add a cache variabie to remove dependency to qpOASES
SET(USE_QPOASES TRUE CACHE BOOL "Use qpOASES solver for static stability")
IF (USE_QPOASES)
  SET(${PROJECT_NAME}_HEADERS ${${PROJECT_NAME}_HEADERS}
    include/hpp/constraints/static-stability.hh
    include/hpp/constraints/qp-static-stability.hh
    )
  PKG_CONFIG_APPEND_CFLAGS (-DHPP_CONSTRAINTS_USE_QPOASES)
ENDIF ()

SETUP_PROJECT()
SET(Eigen_REQUIRED "eigen3 >= 3.2.4")
SEARCH_FOR_EIGEN()
ADD_REQUIRED_DEPENDENCY("hpp-pinocchio >= 4")
ADD_REQUIRED_DEPENDENCY("hpp-statistics")
IF (${USE_QPOASES})
  ADD_REQUIRED_DEPENDENCY("qpOASES >= 3.2")
ENDIF ()

ADD_SUBDIRECTORY (src)

IF (RUN_TESTS)
  SET(BOOST_COMPONENT math unit_test_framework)
  SEARCH_FOR_BOOST()
  ADD_SUBDIRECTORY(tests)
ENDIF ()

PKG_CONFIG_APPEND_LIBS("hpp-constraints")

SETUP_PROJECT_FINALIZE()
