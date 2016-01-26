// Copyright (C) 2016 by Benjamin Chrétien, CNRS-AIST JRL.
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

#ifndef ROBOPTIM_SHARED_TESTS_UTIL_HH
# define ROBOPTIM_SHARED_TESTS_UTIL_HH

// Serialization with Boost
# include <boost/archive/text_oarchive.hpp>
# include <boost/archive/text_iarchive.hpp>

# include "serialize.hh"

typedef boost::filesystem::path path_t;

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

# endif //! ROBOPTIM_SHARED_TESTS_UTIL_HH
