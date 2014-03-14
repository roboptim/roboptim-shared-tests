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

#ifndef ROBOPTIM_SHARED_TESTS_COMMON_HH
# define ROBOPTIM_SHARED_TESTS_COMMON_HH
# include <cmath>
# include <iostream>

# include <boost/mpl/vector.hpp>
# include <boost/mpl/push_back.hpp>
# include <boost/make_shared.hpp>
# include <boost/shared_ptr.hpp>
# include <boost/test/test_case_template.hpp>
# include <boost/test/unit_test.hpp>

# include <log4cxx/xml/domconfigurator.h>

# include <roboptim/core/twice-differentiable-function.hh>
# include <roboptim/core/io.hh>
# include <roboptim/core/numeric-linear-function.hh>
# include <roboptim/core/linear-function.hh>
# include <roboptim/core/optimization-logger.hh>
# include <roboptim/core/solver.hh>
# include <roboptim/core/solver-factory.hh>

# ifndef SOLVER_NAME
#  error "please define solver name"
# endif //! PROBLEM_TYPE

# ifndef PLUGIN_PATH
#  error "please define plug-in path"
# endif //! PROBLEM_TYPE

# ifndef FUNCTION_TYPE
#  error "please define function type"
# endif //! PROBLEM_TYPE

# ifndef TESTS_DATA_DIR
#  error "please define TESTS_DATA_DIR"
# endif //! PROBLEM_TYPE

# ifndef COST_FUNCTION_TYPE
#  error "please define COST_FUNCTION_TYPE"
# endif //! PROBLEM_TYPE

# ifdef LOG_FILENAME
#  define SET_LOG_FILE(solver)						\
  solver.parameters()["ipopt.output_file"].value = std::string(LOG_FILENAME); \
  solver.parameters()["cmaes.output_file"].value = std::string(LOG_FILENAME)
# else //! LOG_FILENAME
#  define SET_LOG_FILE(solver)
# endif //! LOG_FILENAME

typedef FUNCTION_TYPE functionType_t;

// Build extensible constraint type.
typedef boost::mpl::vector< > constraints0_t;

# ifdef CONSTRAINT_TYPE_1
typedef boost::mpl::push_back<
  constraints0_t,
  CONSTRAINT_TYPE_1<functionType_t> >::type
constraints1_t;
# else
typedef constraints0_t constraints1_t;
# endif // CONSTRAINT_TYPE_1

# ifdef CONSTRAINT_TYPE_2
typedef boost::mpl::push_back<
  constraints1_t,
  CONSTRAINT_TYPE_2<functionType_t> >::type
constraints2_t;
# else
typedef constraints1_t constraints2_t;
# endif // CONSTRAINT_TYPE_2

typedef constraints2_t constraints_t;


// Define solver type.
typedef ::roboptim::Solver<COST_FUNCTION_TYPE<functionType_t>, constraints_t >
solver_t;

// See: http://stackoverflow.com/a/20050381/1043187
#define BOOST_CHECK_SMALL_OR_CLOSE(expected, observed, tol)	\
  if (std::fabs (expected) < tol) {				\
    BOOST_CHECK_SMALL(observed, tol);				\
  } else {							\
    BOOST_CHECK_CLOSE(expected, observed, tol);			\
  }

#define CHECK_RESULT(RESULT_TYPE)					\
  /* Get the result. */							\
  RESULT_TYPE& result = boost::get<RESULT_TYPE> (res);			\
  /* Check final x. */							\
  for (F<functionType_t>::size_type i = 0; i < result.x.size (); ++i)	\
    BOOST_CHECK_SMALL_OR_CLOSE (result.x[i], ExpectedResult::x[i], x_tol); \
  /* Check final value. */						\
  BOOST_CHECK_SMALL_OR_CLOSE (result.value[0], ExpectedResult::fx, f_tol); \
  /* Display the result. */						\
  std::cout << "A solution has been found: " << std::endl		\
  << result << std::endl;

#define PROCESS_RESULT()						\
  /* Process the result */						\
  switch (res.which ())							\
    {									\
    case solver_t::SOLVER_VALUE:					\
      {									\
	CHECK_RESULT (Result);						\
	break;								\
      }									\
    case solver_t::SOLVER_VALUE_WARNINGS:				\
      {									\
	CHECK_RESULT (ResultWithWarnings);				\
	break;								\
      }									\
    case solver_t::SOLVER_NO_SOLUTION:					\
      {									\
	std::cout << "A solution should have been found. Failing..."	\
		  << std::endl						\
		  << "No solution was found."				\
		  << std::endl;						\
	BOOST_CHECK_EQUAL (res.which (), solver_t::SOLVER_VALUE);	\
	return;								\
      }									\
    case solver_t::SOLVER_ERROR:					\
      {									\
	std::cout << "A solution should have been found. Failing..."	\
		  << std::endl						\
		  << boost::get<SolverError> (res).what ()		\
		  << std::endl;						\
	BOOST_CHECK_EQUAL (res.which (), solver_t::SOLVER_VALUE);	\
	return;								\
      }									\
    }

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
    lt_dlexit ();
  }
};

#endif //! ROBOPTIM_SHARED_TESTS_COMMON_HH
