# Copyright 2015, Benjamin Chrétien, CNRS-LIRMM.
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

# Benchmark problems to compile.
IF(NOT DEFINED BENCHMARK_PROBLEMS)
  SET(BENCHMARK_PROBLEMS benchmark_17)
ENDIF()

# Benchmark problems that are expected to fail with a given solver.
IF(NOT DEFINED BENCHMARK_PROBLEMS_FAIL)
  #SET(BENCHMARK_PROBLEMS_FAIL)
ENDIF()

# Build benchmark problems.
FOREACH(PROBLEM ${BENCHMARK_PROBLEMS})
  BUILD_TEST("benchmark/${PROBLEM}")
ENDFOREACH()

IF(DEFINED BENCHMARK_PROBLEMS_FAIL)
  MESSAGE("-- WARNING: expecting failure for the following benchmark tests:")
  MESSAGE("            ${BENCHMARK_PROBLEMS_FAIL}")
ENDIF()

# Set WILL_FAIL property for problems that are known to fail.
FOREACH(PROBLEM ${BENCHMARK_PROBLEMS_FAIL})
  # Check that the test is in the list of problems to be compiled.
  LIST(FIND BENCHMARK_PROBLEMS ${PROBLEM} ValidTest)
  IF(NOT ${ValidTest} EQUAL -1)
    EXPECT_TEST_FAIL("benchmark/${PROBLEM}")
  ENDIF()
ENDFOREACH()
