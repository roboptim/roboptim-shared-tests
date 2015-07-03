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
    namespace problem56
    {
      const double a = std::asin (std::sqrt (1./4.2));
      const double b = std::asin (std::sqrt (5./7.2));
      const double c = std::asin (std::sqrt (4./7.));
      const double d = std::asin (std::sqrt (2./7.));

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
	  (7, 1, "-x₀x₁x₂")
      {}

      template <typename T>
      void
      F<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result[0] = -x[0] * x[1] * x[2];
      }

      template <>
      void
      F<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.coeffRef (0) = -x[1] * x[2];
	grad.coeffRef (1) = -x[0] * x[2];
	grad.coeffRef (3) = -x[0] * x[1];
      }

      template <typename T>
      void
      F<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.setZero ();

	grad (0) = -x[1] * x[2];
	grad (1) = -x[0] * x[2];
	grad (3) = -x[0] * x[1];
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
	impl_gradient (gradient_ref, const_argument_ref, size_type)
	  const {}
	void
	impl_jacobian (jacobian_ref jac, const_argument_ref x)
	  const;
      };

      template <typename T>
      G<T>::G ()
	: GenericDifferentiableFunction<T>
	  (7, 4, "x₀ - 4.2sin²x₃, x₁ - 4.2sin²x₄, x₂ - 4.2sin²x₅, "
	   "x₀ + 2x₁ + 2x₂ - 7.2sin²x₆")
      {}

      template <typename T>
      void
      G<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result.setZero ();

	result[0] = x[0] - 4.2 * std::sin (x[3]) * std::sin (x[3]);
	result[1] = x[1] - 4.2 * std::sin (x[4]) * std::sin (x[4]);
	result[2] = x[2] - 4.2 * std::sin (x[5]) * std::sin (x[5]);
	result[3] = x[0] + 2 * x[1] + 2 * x[2]
	  - 7.2 * std::sin (x[6]) * std::sin (x[6]);
      }

      template <>
      void
      G<EigenMatrixSparse>::impl_jacobian
      (jacobian_ref jac, const_argument_ref x) const
      {
	jac.coeffRef (0,0) = 1;
	jac.coeffRef (0,3) = -8.4 * std::sin (x[3]) * std::cos (x[3]);

	jac.coeffRef (1,1) = 1;
	jac.coeffRef (1,4) = -8.4 * std::sin (x[4]) * std::cos (x[4]);

	jac.coeffRef (2,2) = 1;
	jac.coeffRef (2,5) = -8.4 * std::sin (x[5]) * std::cos (x[5]);

	jac.coeffRef (3,0) = 1;
	jac.coeffRef (3,1) = 2;
	jac.coeffRef (3,2) = 2;
	jac.coeffRef (3,6) = -14.4 * std::sin (x[6]) * std::cos (x[6]);
      }

      template <typename T>
      void
      G<T>::impl_jacobian
      (jacobian_ref jac, const_argument_ref x) const
      {
	jac.setZero ();

	jac (0,0) = 1;
	jac (0,3) = -8.4 * std::sin (x[3]) * std::cos (x[3]);

	jac (1,1) = 1;
	jac (1,4) = -8.4 * std::sin (x[4]) * std::cos (x[4]);

	jac (2,2) = 1;
	jac (2,5) = -8.4 * std::sin (x[5]) * std::cos (x[5]);

	jac (3,0) = 1;
	jac (3,1) = 2;
	jac (3,2) = 2;
	jac (3,6) = -14.4 * std::sin (x[6]) * std::cos (x[6]);
      }
    } // end of namespace problem56.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (schittkowski, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (schittkowski_problem56)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem56;

  // Tolerances for Boost checks.
  double f0_tol = 1e-4;
  double x_tol = 1e-4;
  double f_tol = 1e-4;

  ExpectedResult expectedResult;
  expectedResult.f0 = -1.;
  expectedResult.x = (ExpectedResult::argument_t (7)
                      << 2.4, 1.2, 1.2, c, d, d, 0.5 * M_PI).finished ();
  expectedResult.fx = -3.456;

  // Build problem.
  F<functionType_t> f;
  solver_t::problem_t problem (f);

  boost::shared_ptr<G<functionType_t> > g =
    boost::make_shared<G<functionType_t> > ();

  solver_t::problem_t::intervals_t intervals;
  for (F<functionType_t>::size_type i = 0; i < g->outputSize (); ++i)
    intervals.push_back (G<functionType_t>::makeInterval (0., 0.));
  solver_t::problem_t::scaling_t scaling
    (static_cast<std::size_t> (g->outputSize ()), 1.);

  problem.addConstraint (g, intervals, scaling);

  F<functionType_t>::argument_t x (f.inputSize ());
  x << 1, 1, 1, a, a, a, b;
  problem.startingPoint () = x;

  BOOST_CHECK_SMALL_OR_CLOSE (f (x)[0], expectedResult.f0, f0_tol);

  std::cout << f.inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();
  // Set optimization logger
  SET_OPTIMIZATION_LOGGER (solver, "schittkowski/problem-56");

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
