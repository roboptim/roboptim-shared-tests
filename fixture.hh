// Copyright (C) 2013 by Thomas Moulard, AIST, CNRS.
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

#ifndef ROBOPTIM_SHARED_TESTS_FIXTURE_HH
# define ROBOPTIM_SHARED_TESTS_FIXTURE_HH
# include <iostream>
# include <fstream>

# include <boost/make_shared.hpp>
# include <boost/filesystem.hpp>

# include <boost/mpl/list.hpp>
# include <boost/test/output_test_stream.hpp>
# include <boost/test/test_case_template.hpp>
# include <boost/test/unit_test.hpp>

# include <log4cxx/xml/domconfigurator.h>

# include <roboptim/core/numeric-linear-function.hh>
# include <roboptim/core/twice-differentiable-function.hh>
# include <roboptim/core/io.hh>
# include <roboptim/core/solver.hh>
# include <roboptim/core/solver-factory.hh>

struct TestSuiteConfiguration
{
  TestSuiteConfiguration ()
  {
    std::string log4cxxConfigurationFile = TESTS_DATA_DIR;
    log4cxxConfigurationFile += "/log4cxx.xml";
    log4cxx::xml::DOMConfigurator::configure (log4cxxConfigurationFile);

    lt_dlinit();
    BOOST_REQUIRE_EQUAL (lt_dlsetsearchpath (PLUGIN_PATH), 0);
  }

  ~TestSuiteConfiguration ()
  {
    if (lt_dlexit ())
      std::cerr << "lt_dlexit failed" << std::endl;
  }
};

boost::shared_ptr<boost::test_tools::output_test_stream>
retrievePattern (const std::string& testName)
{
  std::string patternFilename = TESTS_DATA_DIR;
  patternFilename += "/";
  patternFilename += testName;
  patternFilename += ".stdout";

  boost::shared_ptr<boost::test_tools::output_test_stream>
    output = boost::make_shared<boost::test_tools::output_test_stream>
    (patternFilename, true);
  return output;
}

#endif //! ROBOPTIM_SHARED_TESTS_FIXTURE_HH
