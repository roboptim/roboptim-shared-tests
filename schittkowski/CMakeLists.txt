# Copyright 2013, Thomas Moulard, CNRS-AIST JRL.
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

IF(NOT DEFINED SCHITTKOWSKI_PROBLEMS)
  # Draft: 11 15 16 20 21 22 24 25 28 33
  SET(SCHITTKOWSKI_PROBLEMS
    1 2 3 4 5 6 7 8 9
    10 12 13 14 17 18 19
    23 26 27 29
    30 31 32 34
    40 41 42 45 47 48 49
    50 51
    71 71b)
ENDIF()

FOREACH(PROBLEM ${SCHITTKOWSKI_PROBLEMS})
  BUILD_TEST("schittkowski/problem_${PROBLEM}")
ENDFOREACH()
