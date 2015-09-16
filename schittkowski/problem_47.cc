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
    namespace problem47
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
	  (5, 1, "(x₀ - x₁)² + (x₁ - x₂)³ + (x₂ - x₃)⁴ + (x₃ - x₄)⁴")
      {}

      template <typename T>
      void
      F<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result[0] = std::pow (x[0] - x[1], 2)
	  + std::pow (x[1] - x[2], 3)
	  + std::pow (x[2] - x[3], 4)
	  + std::pow (x[3] - x[4], 4);
      }

      template <>
      void
      F<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.coeffRef (0) =  2 * (x[0] - x[1]);
	grad.coeffRef (1) = -2 * (x[0] - x[1]) + 3 * std::pow (x[1] - x[2], 2);
	grad.coeffRef (2) = -3 * std::pow (x[1] - x[2], 2) + 4 * std::pow (x[2] - x[3], 3);
	grad.coeffRef (3) = -4 * std::pow (x[2] - x[3], 3) + 4 * std::pow (x[3] - x[4], 3);
	grad.coeffRef (4) = -4 * std::pow (x[3] - x[4], 3);
      }

      template <typename T>
      void
      F<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad[0] =  2 * (x[0] - x[1]);
	grad[1] = -2 * (x[0] - x[1]) + 3 * std::pow (x[1] - x[2], 2);
	grad[2] = -3 * std::pow (x[1] - x[2], 2) + 4 * std::pow (x[2] - x[3], 3);
	grad[3] = -4 * std::pow (x[2] - x[3], 3) + 4 * std::pow (x[3] - x[4], 3);
	grad[4] = -4 * std::pow (x[3] - x[4], 3);
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
	  (5, 3, "x₀ + x₁² + x₂³, x₁ - x₂² + x₃, x₀x₄")
      {}

      template <typename T>
      void
      G<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result[0] = x[0] + std::pow (x[1], 2) + std::pow (x[2], 3);
	result[1] = x[1] - std::pow (x[2], 2) + x[3];
	result[2] = x[0] * x[4];
      }

      template <>
      void
      G<EigenMatrixSparse>::impl_jacobian
      (jacobian_ref jac, const_argument_ref x) const
      {
	jac.coeffRef (0,0) = 1;
	jac.coeffRef (0,1) = 2 * x[1];
	jac.coeffRef (0,2) = 2 * x[2];

	jac.coeffRef (1,1) = 1;
	jac.coeffRef (1,2) = -2 * x[2];
	jac.coeffRef (1,3) = 1;

	jac.coeffRef (2,0) = x[4];
	jac.coeffRef (2,4) = x[0];
      }

      template <typename T>
      void
      G<T>::impl_jacobian
      (jacobian_ref jac, const_argument_ref x) const
      {
	jac.setZero ();

	jac (0,0) = 1;
	jac (0,1) = 2 * x[1];
	jac (0,2) = 2 * x[2];

	jac (1,1) = 1;
	jac (1,2) = -2 * x[2];
	jac (1,3) = 1;

	jac (2,0) = x[4];
	jac (2,4) = x[0];
      }
    } // end of namespace problem47.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (schittkowski, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (schittkowski_problem47)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem47;

  // Tolerances for Boost checks.
  double f0_tol = 1e-4;
  double x_tol = 1e-4;
  double f_tol = 1e-4;

  ExpectedResult expectedResult;
  expectedResult.f0 = 20.7380774886;
  expectedResult.x = (ExpectedResult::argument_t (5) << 1, 1, 1, 1, 1).finished ();
  expectedResult.fx = 0;

  // Build problem.
  F<functionType_t> f;
  solver_t::problem_t problem (f);

  boost::shared_ptr<G<functionType_t> > g =
    boost::make_shared<G<functionType_t> > ();

  solver_t::problem_t::intervals_t intervals;
  intervals.push_back (G<functionType_t>::makeInterval (3., 3.));
  intervals.push_back (G<functionType_t>::makeInterval (1., 1.));
  intervals.push_back (G<functionType_t>::makeInterval (1., 1.));
  solver_t::problem_t::scaling_t scaling
    (static_cast<std::size_t> (g->outputSize ()), 1.);

  problem.addConstraint (g, intervals, scaling);

  F<functionType_t>::argument_t x (5);
  x << 2, std::sqrt(2), -1, 2-std::sqrt(2), 0.5;
  problem.startingPoint () = x;

  BOOST_CHECK_SMALL_OR_CLOSE (f (x)[0], expectedResult.f0, f0_tol);

  std::cout << f.inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();
  // Set optimization logger
  SET_OPTIMIZATION_LOGGER (solver, "schittkowski/problem-47");

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
