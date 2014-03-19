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

# Schittkowski problems to compile.
IF(NOT DEFINED SCHITTKOWSKI_PROBLEMS)
  SET(SCHITTKOWSKI_PROBLEMS
        1  2  3  4  5  6  7  8  9 # Missing problems in PDF:
    10 11 12 13 14 15 16 17 18 19 # ------------------------
    20 21 22 23 24 25 26 27 28 29 #
    30 31 32 33 34 35 36 37 38 39 #
    40 41 42 43    45    47 48 49 # 46
    50 51 52 53    55             # 58
    60 61    63    65 66          #
    71                76          #

    # Particular variants of some problems.
    71b)
ENDIF()

# Schittkowski problems that are expected to fail with a given solver.
IF(NOT DEFINED SCHITTKOWSKI_PROBLEMS_FAIL)
  SET(SCHITTKOWSKI_PROBLEMS_FAIL)
ENDIF()

# Build Schittkowski problems.
FOREACH(PROBLEM ${SCHITTKOWSKI_PROBLEMS})
  BUILD_TEST("schittkowski/problem_${PROBLEM}")
ENDFOREACH()

MESSAGE("-- WARNING: expecting failure for the following Schittkowski tests:")
MESSAGE("            ${SCHITTKOWSKI_PROBLEMS_FAIL}")

# Set WILL_FAIL property for problems that are known to fail.
FOREACH(PROBLEM ${SCHITTKOWSKI_PROBLEMS_FAIL})
  # Check that the test is in the list of problems to be compiled.
  LIST(FIND SCHITTKOWSKI_PROBLEMS ${PROBLEM} ValidTest)
  IF(NOT ${ValidTest} EQUAL -1)
    EXPECT_TEST_FAIL("schittkowski/problem_${PROBLEM}")
  ENDIF()
ENDFOREACH()
