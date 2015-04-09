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

// Serialization with Boost
# include <boost/archive/text_oarchive.hpp>
# include <boost/archive/text_iarchive.hpp>

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

# include "serialize.hh"

typedef boost::filesystem::path path_t;

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

/// \brief Load matrix from a data file.
///
/// This can be used for comparison at a given threshold, rather than
/// relying on the printed matrix.
///
/// \tparam M Eigen matrix type.
/// \param file file containing matrix data. The given relative path
/// should be relative to the tests data directory.
///
/// \return matrix containing the proper data.
template <typename M>
M readMatrix (const path_t& file)
{
  typedef M matrix_t;

  matrix_t m;

  path_t full_path = path_t (TESTS_DATA_DIR) / file;

  std::ifstream ifs (full_path.c_str ());
  boost::archive::text_iarchive ia (ifs);

  ia >> m;

  return m;
}

/// \brief Write matrix to a data file.
///
/// This can be used for later comparison at a given threshold, rather than
/// relying on the printed matrix.
///
/// \tparam M Eigen matrix type.
/// \param file file containing matrix data. The given relative path
/// should be relative to the tests data directory.
template <typename M>
void writeMatrix (const path_t& file, const M& m)
{
  path_t full_path = path_t (TESTS_DATA_DIR) / file;

  std::ofstream ofs (full_path.c_str ());
  boost::archive::text_oarchive oa (ofs);

  oa << m;
}

#endif //! ROBOPTIM_SHARED_TESTS_FIXTURE_HH
