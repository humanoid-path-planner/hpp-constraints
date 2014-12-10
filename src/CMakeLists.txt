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

SET(LIBRARY_NAME ${PROJECT_NAME})
ADD_LIBRARY(${LIBRARY_NAME}
  SHARED
  tool.cc
  position.cc
  orientation.cc
  relative-com.cc
  relative-orientation.cc
  relative-position.cc
  relative-transformation.cc
  static-stability.cc
  )

PKG_CONFIG_USE_DEPENDENCY(${LIBRARY_NAME} hpp-model)

INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION lib)
