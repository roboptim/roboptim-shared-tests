﻿// Copyright (C) 2014 by Thomas Moulard, AIST, CNRS.
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
    namespace problem39
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
	  (4, 1, "-x₀")
      {}

      template <typename T>
      void
      F<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result[0] = -x[0];
      }

      template <>
      void
      F<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref, size_type)
	const
      {
	grad.coeffRef (0) = -1.;
	grad.coeffRef (1) = 0.;
	grad.coeffRef (2) = 0.;
	grad.coeffRef (3) = 0.;
      }

      template <typename T>
      void
      F<T>::impl_gradient (gradient_ref grad, const_argument_ref, size_type)
	const
      {
	grad[0] = -1.;
	grad[1] = 0.;
	grad[2] = 0.;
	grad[3] = 0.;
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
	  (4, 1, "x₁ - x₀³ - x₂²")
      {}

      template <typename T>
      void
      G<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result[0] = x[1] - std::pow (x[0], 3) - std::pow (x[2], 2);
      }

      template <>
      void
      G<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.coeffRef (0) = -3. * std::pow (x[0], 2);
	grad.coeffRef (1) = 1.;
	grad.coeffRef (2) = -2. * x[2];
	grad.coeffRef (3) = 0.;
      }

      template <typename T>
      void
      G<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad[0] = -3. * std::pow (x[0], 2);
	grad[1] = 1.;
	grad[2] = -2. * x[2];
	grad[3] = 0.;
      }

      template <typename T>
      class G2 : public GenericDifferentiableFunction<T>
      {
      public:
	ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
	(GenericDifferentiableFunction<T>);

	explicit G2 ();
	void
	impl_compute (result_ref result, const_argument_ref x) const;
	void
	impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	  const;
      };

      template <typename T>
      G2<T>::G2 ()
	: GenericDifferentiableFunction<T>
	  (4, 1, "x₀² - x₁ - x₃²")
      {}

      template <typename T>
      void
      G2<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result[0] = std::pow (x[0], 2) - x[1] - std::pow (x[3], 2);
      }

      template <>
      void
      G2<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.coeffRef (0) = 2 * x[0];
	grad.coeffRef (1) = -1.;
	grad.coeffRef (2) = 0.;
	grad.coeffRef (3) = -1. * 2. * x[3];
      }

      template <typename T>
      void
      G2<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad[0] = 2 * x[0];
	grad[1] = -1.;
	grad[2] = 0.;
	grad[3] = -1. * 2. * x[3];
      }
    } // end of namespace problem39.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (schittkowski, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (schittkowski_problem39)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem39;

  // Tolerances for Boost checks.
  double f0_tol = 1e-4;
  double x_tol = 1e-4;
  double f_tol = 1e-4;

  ExpectedResult expectedResult;
  expectedResult.f0 = -2.;
  expectedResult.x = (ExpectedResult::argument_t (4) << 1., 1., 0., 0.).finished ();
  expectedResult.fx = -1.;

  // Build problem.
  boost::shared_ptr<F<functionType_t> > f (new F<functionType_t> ());
  solver_t::problem_t problem (f);

  boost::shared_ptr<G<functionType_t> > g =
    boost::make_shared<G<functionType_t> > ();
  problem.addConstraint (g, G<functionType_t>::makeInterval (0., 0.));
  boost::shared_ptr<G2<functionType_t> > g2 =
    boost::make_shared<G2<functionType_t> > ();
  problem.addConstraint (g2, G2<functionType_t>::makeInterval (0., 0.));

  F<functionType_t>::argument_t x (4);
  x << 2., 2., 2., 2.;
  problem.startingPoint () = x;

  BOOST_CHECK_SMALL_OR_CLOSE ((*f) (x)[0], expectedResult.f0, f0_tol);

  std::cout << f->inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();
  // Set optimization logger
  SET_OPTIMIZATION_LOGGER (solver, "schittkowski/problem-39");

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
