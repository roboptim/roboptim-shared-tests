// Copyright (C) 2009 by Thomas Moulard, AIST, CNRS, INRIA.
//
// This file is part of the roboptim.
//
// roboptim is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// roboptim is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with roboptim.  If not, see <http://www.gnu.org/licenses/>.


#ifndef OPTIMIZATION_TESTS_COMMON_HH
# define OPTIMIZATION_TESTS_COMMON_HH
# include <roboptim/core/sys.hh>
# include <roboptim/core/debug.hh>
# include <cstdlib>
# include <fstream>
# include <iomanip>
# include <iostream>
# include <iostream>
# include <stdexcept>
# include <ltdl.h>

# include "config.h"
# include "shared-tests/local-libdir.hh"

static const int TEST_FAILED = 10;
static const int TEST_SUCCEED = 0;

int run_test ();
void init ();
void pauseAtExit ();

void pauseAtExit ()
{
  lt_dlexit();
  lt_dlexit();
#if defined _WIN32 && ROBOPTIM_INTERACTIVE_TESTSUITE
system("PAUSE");
#endif //! _WIN32 && ROBOPTIM_INTERACTIVE_TESTSUITE
}

void init ()
{
#ifdef CWDEBUG
  roboptim::debug::init();
#endif //! CWDEBUG

  std::cout << std::setprecision (2) << std::fixed;

  lt_dlinit();
  if (lt_dlsetsearchpath (LOCAL_LIBDIR))
    std::cerr << "Failed to set search path." << std::endl;

  atexit (pauseAtExit);
}


# define GENERATE_TEST()                                \
  int                                                   \
  main (int argc, char** argv)                          \
  {                                                     \
    init ();						\
							\
    if (argc == 2                                       \
        && std::string (argv[1]) == "--version")        \
      {                                                 \
        std::cout << PACKAGE_STRING << std::endl;       \
        return 0;                                       \
      }                                                 \
                                                        \
    int status = 0;                                     \
    try                                                 \
      {                                                 \
        status = run_test ();                           \
      }                                                 \
    catch (std::runtime_error& e)                       \
      {                                                 \
        std::cerr << e.what () << std::endl;            \
        return 1;                                       \
      }                                                 \
    catch (...)                                         \
      {                                                 \
        std::cerr << "Unexpected error" << std::endl;   \
        return 2;                                       \
      }                                                 \
    return status;                                      \
  }

#define CHECK_FAILURE(EXCEPTION, CMD)           \
  {                                             \
    bool failed = true;                         \
    try                                         \
      {                                         \
        CMD;                                    \
      }                                         \
    catch (EXCEPTION&)                          \
      {                                         \
        failed = false;                         \
      }                                         \
    catch (...)                                 \
      {}                                        \
    if (failed)                                 \
      return TEST_FAILED;                       \
  }

#endif //! OPTIMIZATION_TESTS_COMMON_HH
