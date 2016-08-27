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

# include <roboptim/core/util.hh>

# include "serialize.hh"

typedef boost::filesystem::path path_t;

namespace
{
  template <typename M, typename A>
  inline M genericReadMatrix (const path_t& file)
  {
    typedef M matrix_t;
    typedef A inputArchive_t;

    matrix_t m;

    path_t full_path = path_t (TESTS_DATA_DIR) / file;

    std::ifstream ifs (full_path.c_str ());
    inputArchive_t ia (ifs);

# if (defined ROBOPTIM_HAS_FENV_H && defined ENABLE_SIGFPE)
      // Disable SIGFPE (implementation relies on subnormal numbers)
      roboptim::detail::DisableFPE d;
# endif //! (defined ROBOPTIM_HAS_FENV_H && defined ENABLE_SIGFPE)

    ia >> m;

    return m;
  }

  template <typename M, typename A>
  inline void genericWriteMatrix (const path_t& file, const M& m)
  {
    typedef A ouputArchive_t;

    path_t full_path = path_t (TESTS_DATA_DIR) / file;

    std::ofstream ofs (full_path.c_str ());
    ouputArchive_t oa (ofs);

    oa << m;
  }
} // end of unnamed namespace

/// ** NOTE **
/// We use text archives since binary archives are not portable according to
/// the Boost documentation.

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
  return genericReadMatrix<M, boost::archive::text_iarchive> (file);
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
  genericWriteMatrix<M, boost::archive::text_oarchive> (file, m);
}

# endif //! ROBOPTIM_SHARED_TESTS_UTIL_HH
