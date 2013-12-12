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

#include "common.hh"

namespace roboptim
{
  namespace schittkowski
  {
    namespace problem2
    {
      struct ExpectedResult
      {
	static const double f0;
	static const double x[];
	static const double fx;
      };
      const double ExpectedResult::f0 = 909.;

      // a = (598/1200)^{1/2}
      // b = 400 * a^3
      // x[] = (2 * a * cos (1/3 * arccos(1/b)), 1.5)
      const double ExpectedResult::x[] = {1.224370748736354, 1.5};
      const double ExpectedResult::fx = 0.0504261879;

      template <typename T>
      class F : public GenericDifferentiableFunction<T>
      {
      public:
	ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
	(GenericDifferentiableFunction<T>);

	explicit F () throw ();
	void
	impl_compute (result_t& result, const argument_t& x) const throw ();
	void
	impl_gradient (gradient_t& grad, const argument_t& x, size_type)
	  const throw ();
      };

      template <typename T>
      F<T>::F () throw ()
	: GenericDifferentiableFunction<T>
	  (2, 1, "100 (x₁ - x₀²)² + (1 - x₀)²")
      {}

      template <typename T>
      void
      F<T>::impl_compute (result_t& result, const argument_t& x)
	const throw ()
      {
	result[0] = 100 * std::pow (x[1] - std::pow (x[0], 2), 2)
	  + std::pow (1 - x[0], 2);
      }

      template <>
      void
      F<EigenMatrixSparse>::impl_gradient
      (gradient_t& grad, const argument_t& x, size_type)
	const throw ()
      {
	grad.insert (0) =
	  -400 * x[0] * (x[1] - std::pow (x[0], 2)) - 2 * (1 - x[0]);
	grad.insert (1) = 200 * (x[1] - std::pow (x[0], 2));
      }

      template <typename T>
      void
      F<T>::impl_gradient (gradient_t& grad, const argument_t& x, size_type)
	const throw ()
      {
	grad[0] = -400 * x[0] * (x[1] - std::pow (x[0], 2)) - 2 * (1 - x[0]);
	grad[1] = 200 * (x[1] - std::pow (x[0], 2));
      }
    } // end of namespace problem2.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (schittkowski, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (schittkowski_problem2)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem2;

  // Build problem.
  F<functionType_t> f;
  solver_t::problem_t problem (f);

  // In problem2, this bound changes.
  problem.argumentBounds ()[1] = F<functionType_t>::makeLowerInterval (1.5);

  F<functionType_t>::argument_t x (2);
  x << -2., 1.;
  problem.startingPoint () = x;

  BOOST_CHECK_CLOSE (f (x)[0], ExpectedResult::f0, 1e-6);

  std::cout << f.inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();
  OptimizationLogger<solver_t> logger
    (solver,
     "/tmp/roboptim-shared-tests/" SOLVER_NAME "/schittkowski/problem-2");

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

  // Check if the minimization has succeed.
  if (res.which () != solver_t::SOLVER_VALUE)
    {
      std::cout << "A solution should have been found. Failing..."
                << std::endl
                << boost::get<SolverError> (res).what ()
                << std::endl;
      BOOST_CHECK_EQUAL (res.which (), solver_t::SOLVER_VALUE);
      return;
    }

  // Get the result.
  Result& result = boost::get<Result> (res);

  // Check final x.
  for (F<functionType_t>::vector_t::Index i = 0; i < result.x.size (); ++i)
    BOOST_CHECK_CLOSE (result.x[i], ExpectedResult::x[i], 1e-5);

  // Check final value.
  BOOST_CHECK_CLOSE (1. + result.value[0], 1. + ExpectedResult::fx, 1e-4);


  // Display the result.
  std::cout << "A solution has been found: " << std::endl
	    << result << std::endl;
}

BOOST_AUTO_TEST_SUITE_END ()
