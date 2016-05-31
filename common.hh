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
# ifdef WIN32
#  define _USE_MATH_DEFINES
# endif //WIN32
# include <cmath>
# include <iostream>

# include <boost/ref.hpp>
# include <boost/make_shared.hpp>
# include <boost/mpl/vector.hpp>
# include <boost/mpl/push_back.hpp>
# include <boost/make_shared.hpp>
# include <boost/shared_ptr.hpp>
# include <boost/test/test_case_template.hpp>
# include <boost/test/unit_test.hpp>
# include <boost/test/floating_point_comparison.hpp>

# include <roboptim/core/twice-differentiable-function.hh>
# include <roboptim/core/io.hh>
# include <roboptim/core/numeric-linear-function.hh>
# include <roboptim/core/linear-function.hh>
# include <roboptim/core/optimization-logger.hh>
# include <roboptim/core/solver.hh>
# include <roboptim/core/solver-factory.hh>

# include "fixture.hh"

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

# ifdef LOG_FILENAME
#  define SET_LOG_FILE(solver)						\
  solver.parameters()["ipopt.output_file"].value = std::string(LOG_FILENAME); \
  solver.parameters()["cmaes.output_file"].value = std::string(LOG_FILENAME); \
  solver.parameters()["pagmo.output_file"].value = std::string(LOG_FILENAME); \
  solver.parameters()["knitro.outdir"].value = std::string(LOG_DIR)
# else //! LOG_FILENAME
#  define SET_LOG_FILE(solver)
# endif //! LOG_FILENAME

typedef FUNCTION_TYPE functionType_t;

// Define solver type.
typedef ::roboptim::Solver<functionType_t> solver_t;

typedef ::roboptim::OptimizationLogger<solver_t> logger_t;

namespace roboptim
{
  boost::shared_ptr<logger_t> logger;
} // end of namespace roboptim

#define SET_OPTIMIZATION_LOGGER(SOLVER,FILENAME)	\
  logger = boost::make_shared<logger_t>			\
    (boost::ref<solver_t> (SOLVER),			\
     "/tmp/roboptim-shared-tests/" SOLVER_NAME		\
     "/" FILENAME);

#define RELEASE_OPTIMIZATION_LOGGER()		\
  if (logger)					\
    {						\
      logger.reset ();				\
    }

// Note: tolerances here are in percent, since this is what
// Boost is expecting.

// See: http://stackoverflow.com/a/20050381/1043187
#define BOOST_CHECK_SMALL_OR_CLOSE(EXP, OBS, TOL)	\
  if (std::fabs (EXP) < TOL) {				\
    BOOST_CHECK_SMALL(OBS, TOL);			\
  } else {						\
    BOOST_CHECK_CLOSE(EXP, OBS, TOL);			\
  }

// Run BOOST_SMALL_OR_CLOSE and get result
#define BOOST_SMALL_OR_CLOSE_RES(EXP, OBS, TOL, RES)		\
  if (std::fabs (EXP) < TOL) {					\
    BOOST_CHECK_SMALL (OBS, TOL);				\
    RES = boost::test_tools::check_is_small (OBS, TOL);		\
  } else {							\
    BOOST_CHECK_CLOSE(EXP, OBS, TOL);				\
    RES = (std::fabs (EXP - OBS) < TOL/100. * std::fabs (EXP) \
        && std::fabs (EXP - OBS) < TOL/100. * std::fabs (OBS));	\
  }

// Run BOOST_CHECK and get result
#define BOOST_CHECK_RES(COND, RES)		\
  BOOST_CHECK(COND);				\
  RES = (COND);

// Check the result of the optimization process for an unconstrained problem
#define CHECK_RESULT_UNCONSTRAINED(RESULT_TYPE)				\
  /* Get the result. */							\
  RESULT_TYPE& result = boost::get<RESULT_TYPE> (res);			\
  /* Check final value. */						\
  bool correct_fx = true;						\
  BOOST_SMALL_OR_CLOSE_RES (result.value[0], expectedResult.fx, f_tol, correct_fx); \
  /* Check final bounds on x. */					\
  bool correct_bounds = true;						\
  for (GenericFunction<functionType_t>::size_type i = 0; i < result.x.size (); ++i) { \
    bool res = true;							\
    std::size_t ii = static_cast<std::size_t> (i);			\
    /* Check lower bound. */						\
    BOOST_CHECK_RES (result.x[i]					\
                     - problem.argumentBounds ()[ii].first > -x_tol, res); \
    correct_bounds &= res;						\
    /* Check upper bound. */						\
    BOOST_CHECK_RES (result.x[i]					\
                     - problem.argumentBounds ()[ii].second < x_tol, res); \
    correct_bounds &= res;						\
  }									\
  /* Only check x is we have not found an optimal result. */		\
  if (!(correct_fx && correct_bounds)) {				\
    /* Check final x. */						\
    for (GenericFunction<functionType_t>::size_type i = 0; i < result.x.size (); ++i) \
      BOOST_CHECK_SMALL_OR_CLOSE (result.x[i], expectedResult.x[i], x_tol); \
  }									\
  /* Display the result. */						\
  std::cout << "A solution has been found: " << std::endl		\
  << result << std::endl;


// Check the result of the optimization process
#define CHECK_RESULT(RESULT_TYPE)					\
  /* Get the result. */							\
  RESULT_TYPE& result = boost::get<RESULT_TYPE> (res);			\
  /* Check final value. */						\
  bool success = false;							\
  bool correct_fx = true;						\
  BOOST_SMALL_OR_CLOSE_RES (result.value[0], expectedResult.fx, f_tol, correct_fx); \
  /* Check final bounds on x. */					\
  bool correct_bounds = true;						\
  for (F<functionType_t>::size_type i = 0; i < result.x.size (); ++i) {	\
    bool res = true;							\
    std::size_t ii = static_cast<std::size_t> (i);			\
    /* Check lower bound. */						\
    BOOST_CHECK_RES (result.x[i]					\
                     - problem.argumentBounds ()[ii].first > -x_tol, res); \
    correct_bounds &= res;						\
    /* Check upper bound. */						\
    BOOST_CHECK_RES (result.x[i]					\
                     - problem.argumentBounds ()[ii].second < x_tol, res); \
    correct_bounds &= res;						\
  }									\
  /* Check final constraints. */					\
  bool correct_g = true;						\
  /* Check that final constraint values have been copied to Result. */	\
  F<functionType_t>::size_type n_cstr = 0;				\
  for (size_t i = 0; i < problem.boundsVector ().size (); ++i) {	\
    n_cstr += static_cast<F<functionType_t>::size_type>			\
      (problem.boundsVector ()[i].size ());				\
  }									\
  BOOST_CHECK(n_cstr == static_cast<F<functionType_t>::size_type>	\
              (result.constraints.size ()));				\
  /* For each multidimensional constraint. */				\
  Function::vector_t::Index cstr_i = 0;					\
  for (size_t i = 0; i < problem.boundsVector ().size (); ++i) {	\
    /* For each dimension of the constraint. */				\
    for (size_t j = 0; j < problem.boundsVector ()[i].size (); ++j) {	\
      bool res = true;							\
      /* Check lower bound. */						\
      BOOST_CHECK_RES (result.constraints[cstr_i]			\
                       - problem.boundsVector ()[i][j].first > -f_tol, res); \
      correct_g &= res;							\
      /* Check upper bound. */						\
      BOOST_CHECK_RES (result.constraints[cstr_i]			\
                       - problem.boundsVector ()[i][j].second < f_tol, res); \
      correct_g &= res;							\
      ++cstr_i;								\
    }									\
  }									\
  /* Only check x is we have not found an optimal result. */		\
  if (!(correct_fx && correct_g && correct_bounds)) {			\
    success = true;							\
    /* Check final x. */						\
    for (F<functionType_t>::size_type i = 0; i < result.x.size (); ++i)	\
      {									\
	bool res = false;						\
	BOOST_SMALL_OR_CLOSE_RES (result.x[i], expectedResult.x[i],	\
				  x_tol, res);				\
        success &= res;							\
      }									\
  }									\
  else success = true;							\
  if (logger)								\
    {									\
      if (success)							\
	(*logger) << log_result_true;					\
      else								\
	(*logger) << log_result_false;					\
      (*logger) << solver;						\
    }									\
  /* Display the result. */						\
  std::cout << "A solution has been found: " << std::endl		\
  << result << std::endl;

// Process the result for a constrained problem
#define PROCESS_RESULT()						\
  std::string log_result_true  = "Optimal solution found: true";	\
  std::string log_result_false = "Optimal solution found: false";	\
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
	if (logger)							\
	  {								\
	    (*logger) << log_result_false				\
	              << solver;					\
	    logger.reset ();						\
	  }								\
	return;								\
      }									\
    case solver_t::SOLVER_ERROR:					\
      {									\
	std::cout << "A solution should have been found. Failing..."	\
		  << std::endl						\
		  << boost::get<SolverError> (res).what ()		\
		  << std::endl;						\
	BOOST_CHECK_EQUAL (res.which (), solver_t::SOLVER_VALUE);	\
	if (logger)							\
	  {								\
	    (*logger) << log_result_false				\
	              << solver;					\
	    logger.reset ();						\
	  }								\
	return;								\
      }									\
    }									\
  RELEASE_OPTIMIZATION_LOGGER ();

// Process the result for an unconstrained problem
#define PROCESS_RESULT_UNCONSTRAINED()					\
  std::string log_result_true  = "Optimal solution found: true";	\
  std::string log_result_false = "Optimal solution found: false";	\
  /* Process the result */						\
  switch (res.which ())							\
    {									\
    case solver_t::SOLVER_VALUE:					\
      {									\
	CHECK_RESULT_UNCONSTRAINED (Result);				\
	break;								\
      }									\
    case solver_t::SOLVER_VALUE_WARNINGS:				\
      {									\
	CHECK_RESULT_UNCONSTRAINED (ResultWithWarnings);		\
	break;								\
      }									\
    case solver_t::SOLVER_NO_SOLUTION:					\
      {									\
	std::cout << "A solution should have been found. Failing..."	\
		  << std::endl						\
		  << "No solution was found."				\
		  << std::endl;						\
	BOOST_CHECK_EQUAL (res.which (), solver_t::SOLVER_VALUE);	\
	if (logger)							\
	  {								\
	    (*logger) << log_result_false				\
	              << solver;					\
	    logger.reset ();						\
	  }								\
	return;								\
      }									\
    case solver_t::SOLVER_ERROR:					\
      {									\
	std::cout << "A solution should have been found. Failing..."	\
		  << std::endl						\
		  << boost::get<SolverError> (res).what ()		\
		  << std::endl;						\
	BOOST_CHECK_EQUAL (res.which (), solver_t::SOLVER_VALUE);	\
	if (logger)							\
	  {								\
	    (*logger) << log_result_false				\
	              << solver;					\
	    logger.reset ();						\
	  }								\
	return;								\
      }									\
    }									\
  RELEASE_OPTIMIZATION_LOGGER ();

namespace roboptim
{
  struct ExpectedResult
  {
    typedef solver_t::problem_t::function_t::argument_t argument_t;

    double f0;
    argument_t x;
    double fx;
  };
} // end of namespace roboptim

#endif //! ROBOPTIM_SHARED_TESTS_COMMON_HH
