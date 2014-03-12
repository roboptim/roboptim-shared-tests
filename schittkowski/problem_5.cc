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
    namespace problem5
    {
      struct ExpectedResult
      {
	static const double f0;
	static const double x[];
	static const double fx;
      };
      const double ExpectedResult::f0 = 1.;
      const double ExpectedResult::x[] = {-M_PI / 3. + .5, -M_PI / 3. - .5};
      const double ExpectedResult::fx = -.5 * std::sqrt (3) - (M_PI / 3.);

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
	  (2, 1, "sin (x₀ + x₁) + (x₀ - x₁)² - 1.5x₀ + 2.5x₁ + 1")
      {}

      template <typename T>
      void
      F<T>::impl_compute (result_t& result, const argument_t& x)
	const throw ()
      {
	result[0] = std::sin (x[0] + x[1])
	  + std::pow(x[0] - x[1], 2) - 1.5 * x[0] + 2.5 * x[1] + 1;
      }

      template <>
      void
      F<EigenMatrixSparse>::impl_gradient
      (gradient_t& grad, const argument_t& x, size_type)
	const throw ()
      {
	grad.insert (0) = 2 * x[0] - 2 * x[1] + std::cos (x[0] + x[1]) - 1.5;
	grad.insert (1) = -2 * x[0] + 2 * x[1] + std::cos (x[0] + x[1]) + 2.5;
      }

      template <typename T>
      void
      F<T>::impl_gradient (gradient_t& grad, const argument_t& x, size_type)
	const throw ()
      {
	grad[0] = 2 * x[0] - 2 * x[1] + std::cos (x[0] + x[1]) - 1.5;
	grad[1] = -2 * x[0] + 2 * x[1] + std::cos (x[0] + x[1]) + 2.5;
      }
    } // end of namespace problem4.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (schittkowski, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (schittkowski_problem4)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem5;

  // Build problem.
  F<functionType_t> f;
  solver_t::problem_t problem (f);

  problem.argumentBounds ()[0] = F<functionType_t>::makeInterval (-1.5, 4.);
  problem.argumentBounds ()[1] = F<functionType_t>::makeInterval (-3., 3.);

  F<functionType_t>::argument_t x (2);
  x << 0., 0.;
  problem.startingPoint () = x;

  BOOST_CHECK_CLOSE (f (x)[0], ExpectedResult::f0, 1e-4);

  std::cout << f.inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();
  OptimizationLogger<solver_t> logger
    (solver,
     "/tmp/roboptim-shared-tests/" SOLVER_NAME "/schittkowski/problem-4");

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
    BOOST_CHECK_CLOSE (1. + result.x[i], 1. + ExpectedResult::x[i], 1e-4);

  // Check final value.
  BOOST_CHECK_CLOSE (1. + result.value[0], 1. + ExpectedResult::fx, 1e-4);


  // Display the result.
  std::cout << "A solution has been found: " << std::endl
	    << result << std::endl;
}

BOOST_AUTO_TEST_SUITE_END ()
