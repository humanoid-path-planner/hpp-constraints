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
SET (${LIBRARY_NAME}_SOURCES
  differentiable-function.cc
  distance-between-bodies.cc
  distance-between-points-in-bodies.cc
  tools.cc
  position.cc
  orientation.cc
  transformation.cc
  relative-com.cc
  com-between-feet.cc
  relative-orientation.cc
  relative-position.cc
  relative-transformation.cc
  generic-transformation.cc
  configuration-constraint.cc
  convex-shape-contact.cc
  static-stability.cc
)

IF (${USE_QPOASES})
SET (${LIBRARY_NAME}_SOURCES
  ${${LIBRARY_NAME}_SOURCES}
  qp-static-stability.cc
  )
ENDIF ()

ADD_LIBRARY(${LIBRARY_NAME} SHARED ${${LIBRARY_NAME}_SOURCES})

PKG_CONFIG_USE_DEPENDENCY(${LIBRARY_NAME} hpp-model)
IF (${USE_QPOASES})
  PKG_CONFIG_USE_DEPENDENCY(${LIBRARY_NAME} qpOASES)
ENDIF ()

INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION lib)
