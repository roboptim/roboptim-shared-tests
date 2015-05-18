// Copyright (C) 2015 by Benjamin Chrétien, CNRS-LIRMM.
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

#include <boost/make_shared.hpp>

namespace roboptim
{
  namespace throwTest
  {
    struct ExpectedResult
    {
      static const double x0[];
      static const double fx0;
      static const double x[];
      static const double fx;
    };
    const double ExpectedResult::x0[] = {1.};
    const double ExpectedResult::fx0  = 1.;
    const double ExpectedResult::x[]  = {0.};
    const double ExpectedResult::fx   = 0.;

    /// Function that throws.
    template <typename T>
    struct ThrowF : public GenericDifferentiableFunction<T>
    {
      ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
      (GenericDifferentiableFunction<T>);

      ThrowF () : GenericDifferentiableFunction<T>
	(1, 1, "Throw at the minimum")
      {}

      ~ThrowF ()
      {}

      void impl_compute (result_ref result, const_argument_ref x) const;
      void impl_gradient (gradient_ref gradient, const_argument_ref x,
                          size_type functionId = 0) const;
    };

    template <typename T>
    void
    ThrowF<T>::impl_compute (result_ref result, const_argument_ref x) const
    {
      result.setZero ();

      // Throw a dummy exception
      if (x.isZero ())
        throw std::runtime_error ("null x");

      result[0] = x[0] * x[0];
    }

    template <>
    void
    ThrowF<EigenMatrixSparse>::impl_gradient
    (gradient_ref grad, const_argument_ref x, size_type)
      const
    {
      grad.coeffRef (0) = 2. * x[0];
    }

    template <typename T>
    void
    ThrowF<T>::impl_gradient
    (gradient_ref grad, const_argument_ref x, size_type)
      const
    {
      grad[0] = 2. * x[0];
    }

  } // namespace throwTest
} // namespace roboptim


BOOST_FIXTURE_TEST_SUITE (common, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (throwTest)
{
  using namespace roboptim;
  using namespace roboptim::throwTest;

  // Tolerances for Boost checks.
  double f0_tol = 1e-6;

  // Build problem.
  ThrowF<functionType_t> f;

  solver_t::problem_t problem (f);

  // Load starting point
  ThrowF<functionType_t>::argument_t x (1);
  x << ExpectedResult::x0[0];
  problem.startingPoint () = x;

  // Bounds on x₀
  problem.argumentBounds ()[0] = Function::makeInterval (-1., 1.);

  BOOST_CHECK_SMALL_OR_CLOSE (f (x)[0], ExpectedResult::fx0, f0_tol);

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();

  // Set optimization logger
  SET_OPTIMIZATION_LOGGER (solver, "common/throw");

  // Set optional log file for debugging
  SET_LOG_FILE (solver);

  // Compute the minimum and catch the exception thrown when the minimum is
  //reached.
  BOOST_CHECK_THROW (solver_t::result_t res = solver.minimum (),
                     std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END ()
