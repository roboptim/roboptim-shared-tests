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

IF(SOLVER_NAME)
ELSE()
  MESSAGE(FATAL_ERROR "require SOLVER_NAME to be set")
ENDIF()

IF(PLUGIN_PATH)
ELSE()
  MESSAGE(FATAL_ERROR "require PLUGIN_PATH to be set")
ENDIF()

IF(FUNCTION_TYPE)
ELSE()
  MESSAGE(FATAL_ERROR "require FUNCTION_TYPE to be set")
ENDIF()

# Add Boost path to include directories.
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

# Make Boost.Test generates the main function in test cases.
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

# BUILD_TEST(NAME)
# ----------------
#
# Define a test named `NAME'.
#
# This macro will create a binary from `NAME.cc', link it against
# Boost and add it to the test suite.
#
MACRO(BUILD_TEST NAME)
  ADD_EXECUTABLE(${NAME}${PROGRAM_SUFFIX}
    shared-tests/schittkowski/${NAME}.cc)

  SET_TARGET_PROPERTIES(${NAME}${PROGRAM_SUFFIX}
    PROPERTIES COMPILE_FLAGS
    "-DSOLVER_NAME='\"${SOLVER_NAME}\"' -DPLUGIN_PATH='\"${PLUGIN_PATH}\"' -DFUNCTION_TYPE=${FUNCTION_TYPE}")

  ADD_TEST(${NAME}${PROGRAM_SUFFIX}
    ${RUNTIME_OUTPUT_DIRECTORY}/${NAME}${PROGRAM_SUFFIX})

  PKG_CONFIG_USE_DEPENDENCY(${NAME}${PROGRAM_SUFFIX} roboptim-core)

  # Link against Boost.
  TARGET_LINK_LIBRARIES(${NAME}${PROGRAM_SUFFIX} ${Boost_LIBRARIES})
ENDMACRO()

#BUILD_TEST (distance-to-sphere)

INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/shared-tests/schittkowski/CMakeLists.txt)
