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

SET(CXX_DISABLE_WERROR TRUE)
INCLUDE(cmake/base.cmake)

SET(PROJECT_NAME hpp-constraints)
SET(PROJECT_DESCRIPTION
  "Definition of basic geometric constraints for motion planning")
SET(PROJECT_URL "https://github.com/laas/hpp-constraints")

SET(${PROJECT_NAME}_HEADERS
  include/hpp/constraints/differentiable-function.hh
  include/hpp/constraints/fwd.hh
  include/hpp/constraints/tools.hh
  include/hpp/constraints/orientation.hh
  include/hpp/constraints/position.hh
  include/hpp/constraints/transformation.hh
  include/hpp/constraints/relative-com.hh
  include/hpp/constraints/com-between-feet.hh
  include/hpp/constraints/relative-orientation.hh
  include/hpp/constraints/relative-position.hh
  include/hpp/constraints/relative-transformation.hh
  include/hpp/constraints/static-stability.hh
)

# Activate hpp-util logging if requested
SET (HPP_DEBUG FALSE CACHE BOOL "trigger hpp-util debug output")
IF (HPP_DEBUG)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHPP_DEBUG")
ENDIF()

SETUP_PROJECT()
ADD_REQUIRED_DEPENDENCY("hpp-model >= 3")
PKG_CONFIG_APPEND_LIBS("hpp-constraints")

ADD_SUBDIRECTORY (src)
ADD_SUBDIRECTORY (tests)

SETUP_PROJECT_FINALIZE()
