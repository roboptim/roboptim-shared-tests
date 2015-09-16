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
    namespace problem71
    {
      template <typename T>
      struct F : public GenericTwiceDifferentiableFunction<T>
      {
	ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
	(GenericTwiceDifferentiableFunction<T>);

	F ()
	  : GenericTwiceDifferentiableFunction<T>
	    (4, 1, "x₀ x₃ (x₀ + x₁ + x₂) + x₂")
	{}

	void
	impl_compute (result_ref result, const_argument_ref x) const
	{
	  result[0] = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
	}

	void
	impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	  const;

	void
	impl_hessian (hessian_ref h, const_argument_ref x, size_type)
	  const;
      };

      template <typename T>
      struct G0 : public GenericTwiceDifferentiableFunction<T>
      {
	ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
	(GenericTwiceDifferentiableFunction<T>);

	G0 ()
	  : GenericTwiceDifferentiableFunction<T> (4, 1, "x₀ x₁ x₂ x₃")
	{
	}

	void
	impl_compute (result_ref res, const_argument_ref x) const
	{
	  res.setZero ();
	  res (0) = x[0] * x[1] * x[2] * x[3];
	}

	void
	impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	  const;

	void
	impl_hessian (hessian_ref h, const_argument_ref x, size_type)
	  const;
      };

      template <typename T>
      struct G1 : public GenericTwiceDifferentiableFunction<T>
      {
	ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
	(GenericTwiceDifferentiableFunction<T>);

	G1 ()
	  : GenericTwiceDifferentiableFunction<T>
	    (4, 1, "x₀² + x₁² + x₂² + x₃²")
	{
	}

	void
	impl_compute (result_ref res, const_argument_ref x) const
	{
	  res.setZero ();
	  res (0) = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
	}

	void
	impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	  const;

	void
	impl_hessian (hessian_ref h, const_argument_ref, size_type)
	  const;
      };

      template <>
      void
      F<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.setZero ();
	grad.coeffRef (0) = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
	grad.coeffRef (1) = x[0] * x[3];
	grad.coeffRef (2) = x[0] * x[3] + 1;
	grad.coeffRef (3) = x[0] * (x[0] + x[1] + x[2]);
      }

      template <typename T>
      void
      F<T>::impl_gradient
      (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.setZero ();
	grad[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
	grad[1] = x[0] * x[3];
	grad[2] = x[0] * x[3] + 1;
	grad[3] = x[0] * (x[0] + x[1] + x[2]);
      }

      template <>
      void
      G0<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.reserve (4);
	grad.coeffRef (0) = x[1] * x[2] * x[3];
	grad.coeffRef (1) = x[0] * x[2] * x[3];
	grad.coeffRef (2) = x[0] * x[1] * x[3];
	grad.coeffRef (3) = x[0] * x[1] * x[2];
      }

      template <typename T>
      void
      G0<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.setZero ();
	grad[0] = x[1] * x[2] * x[3];
	grad[1] = x[0] * x[2] * x[3];
	grad[2] = x[0] * x[1] * x[3];
	grad[3] = x[0] * x[1] * x[2];
      }


      template <>
      void
      G1<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.reserve (4);
	grad.coeffRef (0) = 2 * x[0];
	grad.coeffRef (1) = 2 * x[1];
	grad.coeffRef (2) = 2 * x[2];
	grad.coeffRef (3) = 2 * x[3];
      }

      template <typename T>
      void
      G1<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
	const
      {
	grad.setZero ();
	grad[0] = 2 * x[0];
	grad[1] = 2 * x[1];
	grad[2] = 2 * x[2];
	grad[3] = 2 * x[3];
      }

      template <>
      void
      F<EigenMatrixSparse>::impl_hessian
      (hessian_ref h, const_argument_ref x, size_type) const
      {
	h.setZero ();
	h.coeffRef (0, 0) = 2 * x[3];
	h.coeffRef (0, 1) = x[3];
	h.coeffRef (0, 2) = x[3];
	h.coeffRef (0, 3) = 2 * x[0] + x[1] + x[2];

	h.coeffRef (1, 0) = x[3];
	h.coeffRef (1, 3) = x[0];

	h.coeffRef (2, 0) = x[3];
	h.coeffRef (2, 3) = x[1];

	h.coeffRef (3, 0) = 2 * x[0] + x[1] + x[2];
	h.coeffRef (3, 1) = x[0];
	h.coeffRef (3, 2) = x[0];
      }

      template <typename T>
      void
      F<T>::impl_hessian (hessian_ref h, const_argument_ref x, size_type)
	const
      {
	h.setZero ();
	h (0, 0) = 2 * x[3];
	h (0, 1) = x[3];
	h (0, 2) = x[3];
	h (0, 3) = 2 * x[0] + x[1] + x[2];

	h (1, 0) = x[3];
	h (1, 1) = 0.;
	h (1, 2) = 0.;
	h (1, 3) = x[0];

	h (2, 0) = x[3];
	h (2, 1) = 0.;
	h (2, 2) = 0.;
	h (2, 3) = x[1];

	h (3, 0) = 2 * x[0] + x[1] + x[2];
	h (3, 1) = x[0];
	h (3, 2) = x[0];
	h (3, 3) = 0.;
      }

      template <>
      void
      G0<EigenMatrixSparse>::impl_hessian
      (hessian_ref h, const_argument_ref x, size_type) const
      {
	h.setZero ();
	h.coeffRef (0, 1) = x[2] * x[3];
	h.coeffRef (0, 2) = x[1] * x[3];
	h.coeffRef (0, 3) = x[1] * x[2];

	h.coeffRef (1, 0) = x[2] * x[3];
	h.coeffRef (1, 2) = x[0] * x[3];
	h.coeffRef (1, 3) = x[0] * x[2];

	h.coeffRef (2, 0) = x[1] * x[3];
	h.coeffRef (2, 1) = x[0] * x[3];
	h.coeffRef (2, 3) = x[0] * x[1];

	h.coeffRef (3, 0) = x[1] * x[2];
	h.coeffRef (3, 1) = x[0] * x[2];
	h.coeffRef (3, 2) = x[0] * x[1];
      }


      template <typename T>
      void
      G0<T>::impl_hessian (hessian_ref h, const_argument_ref x, size_type)
	const
      {
	h.setZero ();
	h (0, 0) = 0.;
	h (0, 1) = x[2] * x[3];
	h (0, 2) = x[1] * x[3];
	h (0, 3) = x[1] * x[2];

	h (1, 0) = x[2] * x[3];
	h (1, 1) = 0.;
	h (1, 2) = x[0] * x[3];
	h (1, 3) = x[0] * x[2];

	h (2, 0) = x[1] * x[3];
	h (2, 1) = x[0] * x[3];
	h (2, 2) = 0.;
	h (2, 3) = x[0] * x[1];

	h (3, 0) = x[1] * x[2];
	h (3, 1) = x[0] * x[2];
	h (3, 2) = x[0] * x[1];
	h (3, 3) = 0.;
      }


      template <>
      void
      G1<EigenMatrixSparse>::impl_hessian
      (hessian_ref h, const_argument_ref, size_type) const
      {
	h.setZero ();
	h.coeffRef (0, 0) = 2.;
	h.coeffRef (1, 1) = 2.;
	h.coeffRef (2, 2) = 2.;
	h.coeffRef (3, 3) = 2.;
      }

      template <typename T>
      void
      G1<T>::impl_hessian (hessian_ref h, const_argument_ref, size_type)
	const
      {
	h.setZero ();
	h (0, 0) = 2.;
	h (0, 1) = 0.;
	h (0, 2) = 0.;
	h (0, 3) = 0.;

	h (1, 0) = 0.;
	h (1, 1) = 2.;
	h (1, 2) = 0.;
	h (1, 3) = 0.;

	h (2, 0) = 0.;
	h (2, 1) = 0.;
	h (2, 2) = 2.;
	h (2, 3) = 0.;

	h (3, 0) = 0.;
	h (3, 1) = 0.;
	h (3, 2) = 0.;
	h (3, 3) = 2.;
      }
    } // end of namespace problem71.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (schittkowski, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (problem_71)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem71;

  // Tolerances for Boost checks.
  double f0_tol = 1e-6;
  double x_tol = 1e-3;
  double f_tol = 1e-3;

  ExpectedResult expectedResult;
  expectedResult.f0 = 16.;
  expectedResult.x = (ExpectedResult::argument_t (4)
                      << 1., 4.742994, 3.8211503, 1.3794082).finished ();
  expectedResult.fx = 17.0140173;

  // Build problem.
  F<functionType_t> f;
  solver_t::problem_t problem (f);

  // Set bound for all variables.
  // 1. < x_i < 5. (x_i in [1.;5.])
  for (std::size_t i = 0;
       i < static_cast<std::size_t> (problem.function ().inputSize ()); ++i)
    problem.argumentBounds ()[i] = Function::makeInterval (1., 5.);

  // Add constraints.
  boost::shared_ptr<G0<functionType_t> > g0 (new G0<functionType_t> ());
  boost::shared_ptr<G1<functionType_t> > g1 (new G1<functionType_t> ());

  problem.addConstraint
    (boost::static_pointer_cast<
      GenericDifferentiableFunction<functionType_t>  > (g0),
     Function::makeLowerInterval (25.));
  problem.addConstraint
    (boost::static_pointer_cast<
      GenericDifferentiableFunction<functionType_t> > (g1),
     Function::makeInterval (40., 40.));

  // Set the starting point.
  F<functionType_t>::argument_t x (4);
  x << 1., 5., 5., 1.;
  problem.startingPoint () = x;

  BOOST_CHECK_SMALL_OR_CLOSE (f (x)[0], expectedResult.f0, f0_tol);

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();
  // Set optimization logger
  SET_OPTIMIZATION_LOGGER (solver, "schittkowski/problem-71");

  // Set optional log file for debugging
  SET_LOG_FILE(solver);

  // Compute the minimum and retrieve the result.
  solver_t::result_t res = solver.minimum ();

  // Display solver information.
  std::cout << solver << std::endl;

  // Process the result
  PROCESS_RESULT();
}

BOOST_AUTO_TEST_SUITE_END ()
