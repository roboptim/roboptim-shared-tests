# Copyright 2015, Stanislas Brossette, CNRS.
#
# This file is part of roboptim-core.
# roboptim-core is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# roboptim-core is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Lesser Public License for more details.
# You should have received a copy of the GNU Lesser General Public License
# along with roboptim-core.  If not, see <http://www.gnu.org/licenses/>.

IF(NOT DEFINED MANIFOLD_PROBLEMS)
  SET(MANIFOLD_PROBLEMS GeometricProblem SO3Problem)
ENDIF()

FOREACH(PROBLEM ${MANIFOLD_PROBLEMS})
  BUILD_TEST("manifold/${PROBLEM}")
ENDFOREACH()
