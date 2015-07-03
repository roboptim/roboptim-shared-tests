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

#include <roboptim/core/finite-difference-gradient.hh>

namespace roboptim
{
  namespace schittkowski
  {
    namespace problem54b
    {
      const double r = 0.2;

      const double m0 = 1e4;
      const double m1 = 1;
      const double m2 = 2e6;
      const double m3 = 10;
      const double m4 = 1e-3;
      const double m5 = 1e8;

      const double s0 = 8e3;
      const double s1 = 1;
      const double s2 = 7e6;
      const double s3 = 50;
      const double s4 = 0.05;
      const double s5 = 5e8;

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
	  (6, 1, "-exp(-h(x)/2)")
      {}

      template <typename T>
      void
      F<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	value_type h = 1./(1. - r * r) *
	  (std::pow ((x[0] - m0)/s0, 2)
	   + 2 * r * (x[0] - m0) * (x[1] - m1) / (s0 * s1)
	   + std::pow ((x[1] - m1)/s1, 2))
	  + std::pow ((x[2] - m2)/s2, 2)
	  + std::pow ((x[3] - m3)/s3, 2)
	  + std::pow ((x[4] - m4)/s4, 2)
	  + std::pow ((x[5] - m5)/s5, 2);

	result[0] = -std::exp (-h/2.);
      }

      template <typename T>
      void
      F<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type i)
	const
      {
	//FIXME:
	GenericFiniteDifferenceGradient<T> fd (*this);
	fd.gradient (grad, x, i);
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
	  (6, 1, "((x₀ - µ₀) - 0.2σ₀ + 4000(x₁ - µ₁) - 2000σ₁")
      {}

      template <typename T>
      void
      G<T>::impl_compute (result_ref result, const_argument_ref x)
	const
      {
	result[0] = (x[0] - m0) - 0.2 * s0 + 4e3 * (x[1] - m1) - 2e3 * s1;
      }

      template <>
      void
      G<EigenMatrixSparse>::impl_gradient
      (gradient_ref grad, const_argument_ref, size_type)
	const
      {
	grad.coeffRef (0) = 1.;
	grad.coeffRef (1) = 4e3;
      }

      template <typename T>
      void
      G<T>::impl_gradient (gradient_ref grad, const_argument_ref, size_type)
	const
      {
	grad[0] = 1.;
	grad[1] = 4e3;
	grad[2] = 0.;
	grad[3] = 0.;
	grad[4] = 0.;
	grad[5] = 0.;
      }
    } // end of namespace problem54b.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (schittkowski, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (schittkowski_problem54b)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem54b;

  // Tolerances for Boost checks.
  double f0_tol = 1e-4;
  double x_tol = 1e-4;
  double f_tol = 1e-4;

  expectedResult.f0 = -0.7651;
  expectedResult.x = (ExpectedResult::argument_t (6)
                      << 1.2670e4, 1.2322, 1.9999e6, 10., 1e-3, 1e8
                     ).finished ();
  expectedResult.fx = -0.93676;

  // Build problem.
  F<functionType_t> f;
  solver_t::problem_t problem (f);

  problem.argumentBounds ()[0] = F<functionType_t>::makeInterval (0., 2.e4);
  problem.argumentBounds ()[1] = F<functionType_t>::makeInterval (-10., 10.);
  problem.argumentBounds ()[2] = F<functionType_t>::makeInterval (0., 1.e7);
  problem.argumentBounds ()[3] = F<functionType_t>::makeInterval (0., 20.);
  problem.argumentBounds ()[4] = F<functionType_t>::makeInterval (-1., 1.);
  problem.argumentBounds ()[5] = F<functionType_t>::makeInterval (0., 2.e8);

  boost::shared_ptr<G<functionType_t> > g =
    boost::make_shared<G<functionType_t> > ();
  problem.addConstraint (g, G<functionType_t>::makeInterval (0., 0.));

  F<functionType_t>::argument_t x (f.inputSize ());
  x << 6.e3, 1.5, 4e6, 2., 3e-3, 5e7;
  problem.startingPoint () = x;

  BOOST_CHECK_SMALL_OR_CLOSE (f (x)[0], expectedResult.f0, f0_tol);

  std::cout << f.inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();
  // Set optimization logger
  SET_OPTIMIZATION_LOGGER (solver, "schittkowski/problem-54b");

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
