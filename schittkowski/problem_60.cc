// Copyright (C) 2014 by Benjamin Chretien, CNRS.
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
    namespace problem60
    {
      struct ExpectedResult
      {
	static const double f0;
	static const double x[];
	static const double fx;
      };
      const double ExpectedResult::f0 = 1;
      const double ExpectedResult::x[] = {1.104859024, 1.196674194,
                                          1.535262257};
      const double ExpectedResult::fx = 0.03256820025;

      template <typename T>
      class F : public GenericDifferentiableFunction<T>
      {
      public:
	ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
	(GenericDifferentiableFunction<T>);

	explicit F () throw ();
	void
	impl_compute (result_ref result, const_argument_ref& x) const throw ();
	void
	impl_gradient (gradient_ref grad, const_argument_ref& x, size_type)
	  const throw ();
      };

      template <typename T>
      F<T>::F () throw ()
	: GenericDifferentiableFunction<T>
	  (3, 1, "(x₀ - 1)² + (x₀ - x₁)² + (x₁ - x₂)⁴")
      {}

      template <typename T>
      void
      F<T>::impl_compute (result_ref result, const_argument_ref& x)
	const throw ()
      {
	result[0] = std::pow (x[0] - 1, 2)
	  + std::pow (x[0] - x[1], 2)
	  + std::pow (x[1] - x[2], 4);
      }

      template <>
      void
      F<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref& x, size_type)
	const throw ()
      {
	grad.insert (0) = 2. * (2*x[0] - 1 - x[1]);
	grad.insert (1) = -2. * (x[0] - x[1]) + 4. * std::pow (x[1] - x[2], 3);
	grad.insert (2) = -4. * std::pow (x[1] - x[2], 3);
      }

      template <typename T>
      void
      F<T>::impl_gradient (gradient_ref grad, const_argument_ref& x, size_type)
	const throw ()
      {
	grad (0) = 2. * (2*x[0] - 1 - x[1]);
	grad (1) = -2. * (x[0] - x[1]) + 4. * std::pow (x[1] - x[2], 3);
	grad (2) = -4. * std::pow (x[1] - x[2], 3);
      }

      template <typename T>
      class G : public GenericDifferentiableFunction<T>
      {
      public:
	ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
	(GenericDifferentiableFunction<T>);

	explicit G () throw ();
	void
	impl_compute (result_ref result, const_argument_ref& x) const throw ();
	void
	impl_gradient (gradient_ref, const_argument_ref&, size_type)
	  const throw ();
      };

      template <typename T>
      G<T>::G () throw ()
	: GenericDifferentiableFunction<T>
	  (3, 1, "x₀(1 + x₁²) + x₂⁴")
      {}

      template <typename T>
      void
      G<T>::impl_compute (result_ref result, const_argument_ref& x)
	const throw ()
      {
	result[0] = x[0] * (1 + x[1] * x[1]) + std::pow (x[2], 4);
      }

      template <>
      void
      G<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref& x, size_type)
	const throw ()
      {
	grad.insert (0) = 1 + x[1] * x[1];
	grad.insert (1) = 2 * x[0] * x[1];
	grad.insert (2) = 4 * std::pow (x[2], 3);
      }

      template <typename T>
      void
      G<T>::impl_gradient (gradient_ref grad, const_argument_ref& x, size_type)
	const throw ()
      {
	grad[0] = 1 + x[1] * x[1];
	grad[1] = 2 * x[0] * x[1];
	grad[2] = 4 * std::pow (x[2], 3);
      }

    } // end of namespace problem60.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (schittkowski, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (schittkowski_problem60)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem60;

  // Tolerances for Boost checks.
  double f0_tol = 1e-4;
  double x_tol = 1e-4;
  double f_tol = 1e-4;

  // Build problem.
  F<functionType_t> f;
  solver_t::problem_t problem (f);

  for (F<functionType_t>::size_type i = 0; i < f.inputSize (); ++i)
    problem.argumentBounds ()[static_cast<std::size_t> (i)] =
      F<functionType_t>::makeInterval (-10., 10.);

  boost::shared_ptr<G<functionType_t> > g =
    boost::make_shared<G<functionType_t> > ();

  problem.addConstraint (g, G<functionType_t>::makeInterval
                         (4 + 3 * std::sqrt (2), 4 + 3 * std::sqrt (2)));

  F<functionType_t>::argument_t x (f.inputSize ());
  x << 2, 2, 2;
  problem.startingPoint () = x;

  BOOST_CHECK_SMALL_OR_CLOSE (f (x)[0], ExpectedResult::f0, f0_tol);

  std::cout << f.inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();
  OptimizationLogger<solver_t> logger
    (solver,
     "/tmp/roboptim-shared-tests/" SOLVER_NAME "/schittkowski/problem-60");

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
