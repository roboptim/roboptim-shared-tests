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
    namespace problem55
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
	  (6, 1, "x₀ + 2x₁ + 4x₄ + exp(x₀x₃)")
      {}

      template <typename T>
      void
      F<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result[0] = x[0] + 2 * x[1] + 4 * x[4] + std::exp (x[0] * x[3]);
      }

      template <>
      void
      F<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.coeffRef (0) = 1 + x[3] * std::exp (x[0] * x[3]);
	grad.coeffRef (1) = 2;
	grad.coeffRef (3) = x[0] * std::exp (x[0] * x[3]);
	grad.coeffRef (4) = 4;
      }

      template <typename T>
      void
      F<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.setZero ();
	grad (0) = 1 + x[3] * std::exp (x[0] * x[3]);
	grad (1) = 2;
	grad (3) = x[0] * std::exp (x[0] * x[3]);
	grad (4) = 4;
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
	  (6, 6, "x₀ + 2x₁ + 5x₄ - 6, x₀ + x₁ + x₂ - 3, x₃ + x₄ + x₅ - 2, "
	   "x₀ + x₃ - 1, x₁ + x₄ - 2, x₂ + x₅ - 2")
      {}

      template <typename T>
      void
      G<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result[0] = x[0] + 2 * x[1] + 5 * x[4] - 6;
	result[1] = x[0] + x[1] + x[2] - 3;
	result[2] = x[3] + x[4] + x[5] - 2;
	result[3] = x[0] + x[3] - 1;
	result[4] = x[1] + x[4] - 2;
	result[5] = x[2] + x[5] - 2;
      }

      template <>
      void
      G<EigenMatrixSparse>::impl_jacobian
      (jacobian_ref jac, const_argument_ref) const
      {
	jac.coeffRef (0,0) = 1;
	jac.coeffRef (0,1) = 2;
	jac.coeffRef (0,4) = 5;

	jac.coeffRef (1,0) = 1;
	jac.coeffRef (1,1) = 1;
	jac.coeffRef (1,2) = 1;

	jac.coeffRef (2,3) = 1;
	jac.coeffRef (2,4) = 1;
	jac.coeffRef (2,5) = 1;

	jac.coeffRef (3,0) = 1;
	jac.coeffRef (3,3) = 1;

	jac.coeffRef (4,1) = 1;
	jac.coeffRef (4,4) = 1;

	jac.coeffRef (5,2) = 1;
	jac.coeffRef (5,5) = 1;
      }

      template <typename T>
      void
      G<T>::impl_jacobian
      (jacobian_ref jac, const_argument_ref) const
      {
	jac.setZero ();

	jac (0,0) = 1;
	jac (0,1) = 2;
	jac (0,4) = 5;

	jac (1,0) = 1;
	jac (1,1) = 1;
	jac (1,2) = 1;

	jac (2,3) = 1;
	jac (2,4) = 1;
	jac (2,5) = 1;

	jac (3,0) = 1;
	jac (3,3) = 1;

	jac (4,1) = 1;
	jac (4,4) = 1;

	jac (5,2) = 1;
	jac (5,5) = 1;
      }
    } // end of namespace problem55.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (schittkowski, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (schittkowski_problem55)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem55;

  // Tolerances for Boost checks.
  double f0_tol = 1e-4;
  double x_tol = 1e-4;
  double f_tol = 1e-4;

  ExpectedResult expectedResult;
  expectedResult.f0 = 6.;
  expectedResult.x = (ExpectedResult::argument_t (6) << 0., 4./3., 5./3., 1., 2./3., 1./3.).finished ();
  expectedResult.fx = 19./3.;

  // Build problem.
  F<functionType_t> f;
  solver_t::problem_t problem (f);

  problem.argumentBounds ()[0] = F<functionType_t>::makeInterval (0., 1.);
  problem.argumentBounds ()[1] = F<functionType_t>::makeLowerInterval (0.);
  problem.argumentBounds ()[2] = F<functionType_t>::makeLowerInterval (0.);
  problem.argumentBounds ()[3] = F<functionType_t>::makeInterval (0., 1.);
  problem.argumentBounds ()[4] = F<functionType_t>::makeLowerInterval (0.);
  problem.argumentBounds ()[5] = F<functionType_t>::makeLowerInterval (0.);

  boost::shared_ptr<G<functionType_t> > g =
    boost::make_shared<G<functionType_t> > ();

  solver_t::problem_t::intervals_t intervals;
  for (F<functionType_t>::size_type i = 0; i < g->outputSize (); ++i)
    intervals.push_back (G<functionType_t>::makeInterval (0., 0.));
  solver_t::problem_t::scaling_t scaling
    (static_cast<std::size_t> (g->outputSize ()), 1.);

  problem.addConstraint (g, intervals, scaling);

  F<functionType_t>::argument_t x (f.inputSize ());
  x << 1, 2, 0, 0, 0, 2;
  problem.startingPoint () = x;

  BOOST_CHECK_SMALL_OR_CLOSE (f (x)[0], expectedResult.f0, f0_tol);

  std::cout << f.inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();
  // Set optimization logger
  SET_OPTIMIZATION_LOGGER (solver, "schittkowski/problem-55");

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
