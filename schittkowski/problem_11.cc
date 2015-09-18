
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

// THIS TEST DOES NOT GIVE THE GOOD RESULT, HENCE IT IS DISABLED FOR NOW.

#include "common.hh"

namespace roboptim
{
  namespace schittkowski
  {
    namespace problem11
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
	  (2, 1, "(x₀ - 5)² + x₁² - 25")
      {}

      template <typename T>
      void
      F<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result[0] = (x[0] - 5) * (x[0] - 5) + x[1] * x[1] - 25;
      }

      template <>
      void
      F<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.coeffRef (0) = 2 * x[0] - 10;
	grad.coeffRef (1) = 2 * x[1];
      }

      template <typename T>
      void
      F<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad[0] = 2 * x[0] - 10;
	grad[1] = 2 * x[1];
      }

      template <typename T>
      class G : public GenericDifferentiableFunction<T>
      {
      public:
	ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
	(GenericDifferentiableFunction<T>);

	explicit G ();
	void
	impl_compute (result_ref result, const_argument_ref x) const;
	void
	impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	  const;
      };

      template <typename T>
      G<T>::G ()
	: GenericDifferentiableFunction<T>
	  (2, 1, "-x₀² + x₁")
      {}

      template <typename T>
      void
      G<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result[0] = -x[0] * x[0] + x[1];
      }

      template <>
      void
      G<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.coeffRef (0) = -2 * x[0];
	grad.coeffRef (1) = 1;
      }

      template <typename T>
      void
      G<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad[0] = -2 * x[0];
	grad[1] = 1;
      }
    } // end of namespace problem11.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (schittkowski, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (schittkowski_problem11)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem11;

  // Tolerances for Boost checks.
  double f0_tol = 1e-4;
  double x_tol = 1e-4;
  double f_tol = 1e-4;

  ExpectedResult expectedResult;
  expectedResult.f0 = -24.98;
  double a = 7.5 * std::sqrt (6) + std::sqrt (338.5);
  expectedResult.x = (ExpectedResult::argument_t (2) <<
                      (a - (1. / a)) / std::sqrt (6),
                      ((a * a) - 2 + (1 / (a * a))) / 6
                     ).finished ();
  expectedResult.fx = -8.498464223;

  // Build problem.
  boost::shared_ptr<F<functionType_t> > f (new F<functionType_t> ());
  solver_t::problem_t problem (f);

  boost::shared_ptr<G<functionType_t> > g =
    boost::make_shared<G<functionType_t> > ();
  problem.addConstraint (g, G<functionType_t>::makeLowerInterval (0.));

  F<functionType_t>::argument_t x (2);
  x << 4.9, .1;
  problem.startingPoint () = x;

  BOOST_CHECK_SMALL_OR_CLOSE ((*f) (x)[0], expectedResult.f0, f0_tol);

  std::cout << f->inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();
  // Set optimization logger
  SET_OPTIMIZATION_LOGGER (solver, "schittkowski/problem-11");

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
