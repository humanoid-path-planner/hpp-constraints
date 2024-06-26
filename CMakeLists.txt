# Copyright 2014, 2020, Florent Lamiraux, Guilhem Saurel, CNRS
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

cmake_minimum_required(VERSION 3.10)

set(PROJECT_NAME hpp-constraints)
set(PROJECT_DESCRIPTION
    "Definition of basic geometric constraints for motion planning")

set(CXX_DISABLE_WERROR TRUE)
set(DOXYGEN_USE_MATHJAX YES)
set(PROJECT_USE_CMAKE_EXPORT TRUE)

# Check if the submodule cmake have been initialized
set(JRL_CMAKE_MODULES "${CMAKE_CURRENT_LIST_DIR}/cmake")
if(EXISTS "${JRL_CMAKE_MODULES}/base.cmake")
  message(STATUS "JRL cmakemodules found in 'cmake/' git submodule")
else()
  find_package(jrl-cmakemodules QUIET CONFIG)
  if(jrl-cmakemodules_FOUND)
    get_property(
      JRL_CMAKE_MODULES
      TARGET jrl-cmakemodules::jrl-cmakemodules
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    message(STATUS "JRL cmakemodules found on system at ${JRL_CMAKE_MODULES}")
  elseif(${CMAKE_VERSION} VERSION_LESS "3.14.0")
    message(
      FATAL_ERROR
        "\nCan't find jrl-cmakemodules. Please either:\n"
        "  - use git submodule: 'git submodule update --init'\n"
        "  - or install https://github.com/jrl-umi3218/jrl-cmakemodules\n"
        "  - or upgrade your CMake version to >= 3.14 to allow automatic fetching\n"
    )
  else()
    message(STATUS "JRL cmakemodules not found. Let's fetch it.")
    include(FetchContent)
    FetchContent_Declare(
      "jrl-cmakemodules"
      GIT_REPOSITORY "https://github.com/jrl-umi3218/jrl-cmakemodules.git")
    FetchContent_MakeAvailable("jrl-cmakemodules")
    FetchContent_GetProperties("jrl-cmakemodules" SOURCE_DIR JRL_CMAKE_MODULES)
  endif()
endif()

include("${JRL_CMAKE_MODULES}/hpp.cmake")
include("${JRL_CMAKE_MODULES}/boost.cmake")

# Activate hpp-util logging if requested
set(HPP_CONSTRAINTS_NUMERIC_DEBUG
    ${HPP_DEBUG}
    CACHE BOOL "trigger hpp-util debug output")
if(HPP_CONSTRAINTS_NUMERIC_DEBUG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHPP_CONSTRAINTS_NUMERIC_DEBUG")
endif()

option(CHECK_JACOBIANS "Check jacobians against finite difference." OFF)
if(CHECK_JACOBIANS)
  add_definitions(-DCHECK_JACOBIANS)
endif(CHECK_JACOBIANS)

# Add a cache variable to remove dependency to qpOASES
set(USE_QPOASES
    TRUE
    CACHE BOOL "Use qpOASES solver for static stability")

compute_project_args(PROJECT_ARGS LANGUAGES CXX)
project(${PROJECT_NAME} ${PROJECT_ARGS})

add_project_dependency(hpp-pinocchio REQUIRED)
if(USE_QPOASES)
  set(CMAKE_MODULE_PATH "${JRL_CMAKE_MODULES}/find-external/qpOASES")
  find_package(qpOASES REQUIRED)
endif()

set(${PROJECT_NAME}_HEADERS
    include/hpp/constraints/differentiable-function.hh
    include/hpp/constraints/differentiable-function-set.hh
    include/hpp/constraints/active-set-differentiable-function.hh
    include/hpp/constraints/affine-function.hh
    include/hpp/constraints/comparison-types.hh
    include/hpp/constraints/distance-between-bodies.hh
    include/hpp/constraints/fwd.hh
    include/hpp/constraints/svd.hh
    include/hpp/constraints/tools.hh
    include/hpp/constraints/matrix-view.hh
    include/hpp/constraints/generic-transformation.hh
    include/hpp/constraints/macros.hh
    include/hpp/constraints/manipulability.hh
    include/hpp/constraints/convex-shape.hh
    include/hpp/constraints/convex-shape-contact.hh
    include/hpp/constraints/static-stability.hh
    include/hpp/constraints/symbolic-calculus.hh
    include/hpp/constraints/symbolic-function.hh
    include/hpp/constraints/relative-com.hh
    include/hpp/constraints/com-between-feet.hh
    include/hpp/constraints/configuration-constraint.hh
    include/hpp/constraints/explicit-constraint-set.hh
    include/hpp/constraints/explicit.hh
    include/hpp/constraints/explicit/implicit-function.hh
    include/hpp/constraints/explicit/relative-pose.hh
    include/hpp/constraints/explicit/convex-shape-contact.hh
    include/hpp/constraints/explicit/relative-transformation.hh
    include/hpp/constraints/implicit.hh
    include/hpp/constraints/implicit-constraint-set.hh
    include/hpp/constraints/locked-joint.hh
    include/hpp/constraints/serialization.hh
    include/hpp/constraints/solver/hierarchical-iterative.hh
    include/hpp/constraints/solver/by-substitution.hh
    include/hpp/constraints/function/of-parameter-subset.hh
    include/hpp/constraints/function/difference.hh
    include/hpp/constraints/solver/impl/by-substitution.hh
    include/hpp/constraints/solver/impl/hierarchical-iterative.hh
    include/hpp/constraints/impl/matrix-view.hh
    include/hpp/constraints/impl/matrix-view-operation.hh)

set(${PROJECT_NAME}_SOURCES
    src/affine-function.cc
    src/differentiable-function.cc
    src/differentiable-function-set.cc
    src/generic-transformation.cc
    src/relative-com.cc
    src/com-between-feet.cc
    src/distance-between-bodies.cc
    src/distance-between-points-in-bodies.cc
    src/configuration-constraint.cc
    src/convex-shape-contact.cc
    src/matrix-view.cc
    src/manipulability.cc
    src/static-stability.cc
    src/explicit-constraint-set.cc
    src/implicit.cc
    src/explicit.cc
    src/explicit/convex-shape-contact.cc
    src/explicit/implicit-function.cc
    src/explicit/input-configurations.hh
    src/explicit/relative-transformation.cc
    src/explicit/relative-pose.cc
    src/function/of-parameter-subset.cc
    src/function/difference.cc
    src/locked-joint.cc
    src/solver/by-substitution.cc
    src/solver/hierarchical-iterative.cc)

if(USE_QPOASES)
  set(${PROJECT_NAME}_HEADERS ${${PROJECT_NAME}_HEADERS}
                              include/hpp/constraints/qp-static-stability.hh)
  pkg_config_append_cflags(-DHPP_CONSTRAINTS_USE_QPOASES)
  set(${PROJECT_NAME}_SOURCES ${${PROJECT_NAME}_SOURCES}
                              src/qp-static-stability.cc)
endif(USE_QPOASES)

add_library(${PROJECT_NAME} SHARED ${${PROJECT_NAME}_SOURCES}
                                   ${${PROJECT_NAME}_HEADERS})
target_include_directories(
  ${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(${PROJECT_NAME} PUBLIC hpp-pinocchio::hpp-pinocchio)

if(USE_QPOASES)
  target_include_directories(${PROJECT_NAME} PUBLIC ${qpOASES_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${qpOASES_LIBRARIES})
  target_compile_definitions(${PROJECT_NAME}
                             PUBLIC -DHPP_CONSTRAINTS_USE_QPOASES)
endif(USE_QPOASES)

install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${TARGETS_EXPORT_NAME}
  DESTINATION lib)

if(BUILD_TESTING)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework)
  find_package(example-robot-data REQUIRED)
  add_subdirectory(tests)
endif()

pkg_config_append_libs("hpp-constraints")

install(FILES package.xml DESTINATION share/${PROJECT_NAME})
