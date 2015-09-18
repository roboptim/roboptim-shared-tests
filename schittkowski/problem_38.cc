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
    namespace problem38
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
	  (4, 1,
	   "100 (x₁ - x₀²)² + (1 - x₀)² + 90 (x₃ - x₂²)² + (1 - x₂)²"
	   " + 10.1 ((x₁ - 1)² + (x₃ - 1)²) + 19.8 (x₁ - 1) (x₃ - 1)")
      {}

      template <typename T>
      void
      F<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result[0] =
	  100. * std::pow (x[1] - std::pow (x[0], 2.), 2.)
	  + std::pow (1 - x[0], 2.)
	  + 90. * std::pow (x[3] - std::pow (x[2], 2.), 2.)
	  + std::pow (1. - x[2], 2.)
	  + 10.1 * (std::pow (x[1] - 1., 2.) + std::pow (x[3] - 1., 2.))
	  + 19.8 * (x[1] - 1.) * (x[3] - 1.);
      }

      template <>
      void
      F<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.coeffRef (0) =
	  400. * std::pow (x[0], 3) - 400. * x[0] * x[1] + 2. * x[0] - 2.;
	grad.coeffRef (1) =
	  -200. * std::pow (x[0], 2)  + 220.2 * x[1] + 19.8 * x[3] - 40.;
	grad.coeffRef (2) =
	  360. * std::pow (x[2], 3)  - 360. * x[2] * x[3] + 2. * x[2] - 2.;
	grad.coeffRef (3) =
	  19.8 * x[1] - 180. * std::pow (x[2], 2)  + 200.2 * x[3] - 40.;
      }

      template <typename T>
      void
      F<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad[0] =
	  400. * std::pow (x[0], 3) - 400. * x[0] * x[1] + 2. * x[0] - 2.;
	grad[1] =
	  -200. * std::pow (x[0], 2)  + 220.2 * x[1] + 19.8 * x[3] - 40.;
	grad[2] =
	  360. * std::pow (x[2], 3)  - 360. * x[2] * x[3] + 2. * x[2] - 2.;
	grad[3] =
	  19.8 * x[1] - 180. * std::pow (x[2], 2)  + 200.2 * x[3] - 40.;
      }
    } // end of namespace problem38.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (schittkowski, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (schittkowski_problem38)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem38;

  // Tolerances for Boost checks.
  double f0_tol = 1e-4;
  double x_tol = 1e-4;
  double f_tol = 1e-4;

  ExpectedResult expectedResult;
  expectedResult.f0 = 19192;
  expectedResult.x = (ExpectedResult::argument_t (4) << 1., 1., 1., 1.).finished ();
  expectedResult.fx = 0.;

  // Build problem.
  boost::shared_ptr<F<functionType_t> > f (new F<functionType_t> ());
  solver_t::problem_t problem (f);

  for (std::size_t i = 0; i < 4; ++i)
    problem.argumentBounds ()[i] = F<functionType_t>::makeInterval (-10., 10.);

  F<functionType_t>::argument_t x (4);
  x << -3., -1., -3., -1.;
  problem.startingPoint () = x;

  BOOST_CHECK_SMALL_OR_CLOSE ((*f) (x)[0], expectedResult.f0, f0_tol);

  std::cout << f->inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();
  // Set optimization logger
  SET_OPTIMIZATION_LOGGER (solver, "schittkowski/problem-38");

  // Set optional log file for debugging
  SET_LOG_FILE(solver);

  std::cout << f->inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Compute the minimum and retrieve the result.
  solver_t::result_t res = solver.minimum ();

  std::cout << f->inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Display solver information.
  std::cout << solver << std::endl;

  // Process the result
  PROCESS_RESULT();
}

BOOST_AUTO_TEST_SUITE_END ()
