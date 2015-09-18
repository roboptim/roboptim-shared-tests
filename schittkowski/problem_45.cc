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
    namespace problem45
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
	  (5, 1, "2 - (1/120)x₀x₁x₂x₃x₄")
      {}

      template <typename T>
      void
      F<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result[0] = 2. - (1./120.)*x[0]*x[1]*x[2]*x[3]*x[4];
      }

      template <>
      void
      F<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.coeffRef (0) = -(1./120.)*x[1]*x[2]*x[3]*x[4];
	grad.coeffRef (1) = -(1./120.)*x[0]*x[2]*x[3]*x[4];
	grad.coeffRef (2) = -(1./120.)*x[0]*x[1]*x[3]*x[4];
	grad.coeffRef (3) = -(1./120.)*x[0]*x[1]*x[2]*x[4];
	grad.coeffRef (4) = -(1./120.)*x[0]*x[1]*x[2]*x[3];
      }

      template <typename T>
      void
      F<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad[0] = -(1./120.)*x[1]*x[2]*x[3]*x[4];
	grad[1] = -(1./120.)*x[0]*x[2]*x[3]*x[4];
	grad[2] = -(1./120.)*x[0]*x[1]*x[3]*x[4];
	grad[3] = -(1./120.)*x[0]*x[1]*x[2]*x[4];
	grad[4] = -(1./120.)*x[0]*x[1]*x[2]*x[3];
      }
    } // end of namespace problem45.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (schittkowski, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (schittkowski_problem45)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem45;

  // Tolerances for Boost checks.
  double f0_tol = 1e-4;
  double x_tol = 1e-4;
  double f_tol = 1e-4;

  ExpectedResult expectedResult;
  expectedResult.f0 = 26./15.;
  expectedResult.x = (ExpectedResult::argument_t (5)
                      << 1., 2., 3., 4., 5.).finished ();
  expectedResult.fx = 1.;

  // Build problem.
  boost::shared_ptr<F<functionType_t> > f (new F<functionType_t> ());
  solver_t::problem_t problem (f);

  for (F<functionType_t>::size_type i = 0; i < f->inputSize (); ++i)
    problem.argumentBounds ()[static_cast<std::size_t> (i)]
      = F<functionType_t>::makeInterval (0., 1. + static_cast<double> (i));

  F<functionType_t>::argument_t x (5);
  x << 2, 2, 2, 2, 2;
  problem.startingPoint () = x;

  BOOST_CHECK_SMALL_OR_CLOSE ((*f) (x)[0], expectedResult.f0, f0_tol);

  std::cout << f->inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();
  // Set optimization logger
  SET_OPTIMIZATION_LOGGER (solver, "schittkowski/problem-45");

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
