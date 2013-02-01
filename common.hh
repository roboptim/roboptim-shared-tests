// Copyright (C) 2009-2012 by Thomas Moulard, AIST, CNRS, INRIA.
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
# include <cstdlib>
# include <fstream>
# include <iomanip>
# include <iostream>
# include <iostream>
# include <stdexcept>
# include <string>
# include <ltdl.h>

# include <boost/make_shared.hpp>
# include <boost/shared_ptr.hpp>
# include <boost/test/unit_test.hpp>
# include <boost/test/output_test_stream.hpp> 

# include <log4cxx/basicconfigurator.h>

# include <roboptim/core/util.hh>

# include "config.h"
# include "shared-tests/local-libdir.hh"

void pauseAtExit ();
boost::shared_ptr<boost::test_tools::output_test_stream>
retrievePattern (const std::string& testName);
boost::unit_test::test_suite*
init_unit_test_suite (int argc, char* argv[]);

void pauseAtExit ()
{
  lt_dlexit();
  lt_dlexit();
#if defined _WIN32 && ROBOPTIM_INTERACTIVE_TESTSUITE
system("PAUSE");
#endif //! _WIN32 && ROBOPTIM_INTERACTIVE_TESTSUITE
}

boost::shared_ptr<boost::test_tools::output_test_stream>
retrievePattern (const std::string& testName)
{
  std::string patternFilename = LOCAL_TESTDIR;
  patternFilename += "/";
  patternFilename += testName;
  patternFilename += ".stdout";

  std::cout << patternFilename << std::endl;

  boost::shared_ptr<boost::test_tools::output_test_stream>
    output = boost::make_shared<boost::test_tools::output_test_stream>
    (patternFilename, true);
  return output;
}

boost::unit_test::test_suite*
init_unit_test_suite (int, char*[])
{
  log4cxx::BasicConfigurator::configure ();

  std::cout << std::setprecision (2) << std::fixed;

  lt_dlinit();
  if (lt_dlsetsearchpath (LOCAL_LIBDIR))
    std::cerr << "Failed to set search path." << std::endl;

  atexit (pauseAtExit);
  return 0;
}

#endif //! OPTIMIZATION_TESTS_COMMON_HH
