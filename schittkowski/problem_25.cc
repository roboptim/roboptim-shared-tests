// Copyright (C) 2014 by Thomas Moulard, AIST, CNRS.
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

#include "common.hh"

namespace roboptim
{
  namespace schittkowski
  {
    namespace problem25
    {
      template <typename T>
      class F : public GenericDifferentiableFunction<T>
      {
      public:
	ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
	(GenericDifferentiableFunction<T>);

	explicit F ();
	void
	impl_compute (result_ref result, const_argument_ref x) const;
	void
	impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	  const;
      };

      template <typename T>
      F<T>::F ()
	: GenericDifferentiableFunction<T>
	  (3, 1, "Σ f_i(x)²")
      {}

      template <typename T>
      void
      F<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result[0] = 0.;
	for (int i = 1; i < 100; ++i)
	  {
	    value_type u = 25 + std::pow (-50. * std::log (0.01 * i), 2. / 3.);
	    value_type f =
	      -0.01 * i + std::exp ((- 1. / x[0]) * std::pow (u - x[1], x[2]));
	    result[0] += f * f;
	  }
      }

      template <>
      void
      F<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	value_type dx0 = 0.;
	value_type dx1 = 0.;
	value_type dx2 = 0.;
	for (int i = 1; i < 100; ++i)
	  {
	    value_type u = 25 + std::pow (-50. * std::log (0.01 * i), 2. / 3.);
	    dx0 +=
	      1. / (x[0] * x[0])
	      * std::pow (u - x[1], x[2])
	      * std::exp (-1. / x[0] * std::pow (u - x[1], x[2]));


	    dx1 +=
	      1. / (x[0] * (u - x[1]))
	      * x[2] * std::pow (u - x[1], x[2])
	      * std::exp (-1. / x[0] * std::pow (u - x[1], x[2]));
	    dx2 +=
	      -1. / x[0]
	      * std::pow (u - x[1], x[2])
	      * std::exp (-1. / x[0] * std::pow (u - x[1], x[2]))
	      * std::log (u - x[1]);
	  }


	grad.coeffRef (0) = dx0;
	grad.coeffRef (1) = dx1;
	grad.coeffRef (2) = dx2;
      }

      template <typename T>
      void
      F<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	value_type dx0 = 0.;
	value_type dx1 = 0.;
	value_type dx2 = 0.;
	for (int i = 1; i < 100; ++i)
	  {
	    value_type u = 25 + std::pow (-50. * std::log (0.01 * i), 2. / 3.);
	    dx0 +=
	      1. / (x[0] * x[0])
	      * std::pow (u - x[1], x[2])
	      * std::exp (-1. / x[0] * std::pow (u - x[1], x[2]));


	    dx1 +=
	      1. / (x[0] * (u - x[1]))
	      * x[2] * std::pow (u - x[1], x[2])
	      * std::exp (-1. / x[0] * std::pow (u - x[1], x[2]));
	    dx2 +=
	      -1. / x[0]
	      * std::pow (u - x[1], x[2])
	      * std::exp (-1. / x[0] * std::pow (u - x[1], x[2]))
	      * std::log (u - x[1]);
	  }
	grad[0] = dx0;
	grad[1] = dx1;
  	grad[2] = dx2;
      }
    } // end of namespace problem25.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (schittkowski, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (schittkowski_problem25)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem25;

  // Tolerances for Boost checks.
  double f0_tol = 1e-4;
  double x_tol = 1e-4;
  double f_tol = 1e-4;

  ExpectedResult expectedResult;
  expectedResult.f0 = 32.835;
  expectedResult.x = (ExpectedResult::argument_t (3) << 50., 25., 1.5).finished ();
  expectedResult.fx = 0.;

  // Build problem.
  F<functionType_t> f;
  solver_t::problem_t problem (f);

  problem.argumentBounds ()[0] = F<functionType_t>::makeInterval (0.1, 100.);
  problem.argumentBounds ()[1] = F<functionType_t>::makeInterval (0., 25.6);
  problem.argumentBounds ()[2] = F<functionType_t>::makeInterval (0., 5.);

  F<functionType_t>::argument_t x (3);
  x << 100, 12.5, 3;
  problem.startingPoint () = x;

  BOOST_CHECK_SMALL_OR_CLOSE (f (x)[0], expectedResult.f0, f0_tol);

  std::cout << f.inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();
  // Set optimization logger
  SET_OPTIMIZATION_LOGGER (solver, "schittkowski/problem-25");

  // Set optional log file for debugging
  SET_LOG_FILE(solver);

  std::cout << f.inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Compute the minimum and retrieve the result.
  solver_t::result_t res = solver.minimum ();

  std::cout << f.inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Display solver information.
  std::cout << solver << std::endl;

  // Process the result
  PROCESS_RESULT();
}

BOOST_AUTO_TEST_SUITE_END ()
