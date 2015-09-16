// Copyright (C) 2011-2013 by Florent Lamiraux, Thomas Moulard, AIST, CNRS.
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

#include <roboptim/core/sum-of-c1-squares.hh>

namespace roboptim
{
  namespace distanceToSphere
  {
    /// Distance between a point on unit sphere and another point in R^3
    template <typename T>
    struct F : public GenericDifferentiableFunction<T>
    {
      ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
      (GenericDifferentiableFunction<T>);

      explicit F (const ExpectedResult& target) : GenericDifferentiableFunction<T>
		      (2, 3,
		       "vector between unit sphere and point (x,y,z)"),
		      point_ (3)
      {
	sphericalCoordinates (point_, target.x[0], target.x[1]);
	point_ *= 2.;
      }

      ~F ()
      {}

      void impl_compute (result_ref result, const_argument_ref x) const;
      void impl_gradient (gradient_ref gradient, const_argument_ref x,
                          size_type functionId = 0) const;


      static void sphericalCoordinates (result_ref res, double theta, double phi)
      {
        res (0) = cos(theta) * cos(phi);
        res (1) = sin(theta) * cos(phi);
        res (2) = sin(phi);
      }
      result_t point_;
    };

    template <typename T>
    void
    F<T>::impl_compute (result_ref result, const_argument_ref x) const
    {
      result.setZero ();
      double theta = x[0];
      double phi = x[1];
      sphericalCoordinates (result, theta, phi);
      result -= point_;
    }

    template <>
    void
    F<EigenMatrixSparse>::impl_gradient
    (gradient_ref grad, const_argument_ref x, size_type functionId)
      const
    {
      grad.setZero ();

      double theta = x[0];
      double phi = x[1];

      switch (functionId)
        {
        case 0:
          grad.coeffRef (0) = -sin(theta) * cos(phi);
          grad.coeffRef (1) = -cos(theta) * sin(phi);
          break;
        case 1:
          grad.coeffRef (0) = cos(theta) * cos(phi);
          grad.coeffRef (1) = -sin(theta) * sin(phi);
          break;
        case 2:
          grad.coeffRef (0) = 0.;
          grad.coeffRef (1) = cos(phi);
          break;
        default:
          abort();
        }
    }

    template <typename T>
    void
    F<T>::impl_gradient
    (gradient_ref grad, const_argument_ref x, size_type functionId)
      const
    {
      grad.setZero ();

      double theta = x[0];
      double phi = x[1];

      switch (functionId)
        {
        case 0:
          grad[0] = -sin(theta) * cos(phi);
          grad[1] = -cos(theta) * sin(phi);
          break;
        case 1:
          grad[0] = cos(theta) * cos(phi);
          grad[1] = -sin(theta) * sin(phi);
          break;
        case 2:
          grad[0] = 0.;
          grad[1] = cos(phi);
          break;
        default:
          abort();
        }
    }

  } // namespace distanceToSphere
} // namespace roboptim


BOOST_FIXTURE_TEST_SUITE (distanceToSphere, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (distanceToSphere_problem1)
{
  using namespace roboptim;
  using namespace roboptim::distanceToSphere;

  // Tolerances for Boost checks.
  double f0_tol = 1e-6;
  double x_tol = 1e-5;
  double f_tol = 1e-4;

  ExpectedResult expectedResult;
  expectedResult.f0 = 4.8974713057829096;
  expectedResult.x = (ExpectedResult::argument_t (2) << -1.5, -1.2).finished ();
  expectedResult.fx = 1.0;

  // Build problem.
  boost::shared_ptr <F<functionType_t> > f (new F<functionType_t> (expectedResult));
  GenericSumOfC1Squares<functionType_t> soq (f, "");

  solver_t::problem_t problem (soq);

  // Load starting point
  F<functionType_t>::argument_t x (2);
  x << 0., 0.;
  problem.startingPoint () = x;

  // Set arguments names (optional).
  F<functionType_t>::names_t
    names (static_cast<std::size_t> (f->inputSize ()));
  names[0] = "θ";
  names[1] = "φ";
  problem.argumentNames () = names;

  // Bounds on theta \in [-Pi/2, Pi/2]
  problem.argumentBounds ()[0] = Function::makeInterval (-M_PI_2, M_PI_2);

  // Bounds on phi \in [-Pi, Pi]
  problem.argumentBounds ()[1] = Function::makeInterval (-M_PI, M_PI);

  BOOST_CHECK_SMALL_OR_CLOSE (soq (x)[0], expectedResult.f0, f0_tol);

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();

  // Set optimization logger
  SET_OPTIMIZATION_LOGGER (solver, "roboptim/distance-to-sphere");

  // Set optional log file for debugging
  SET_LOG_FILE (solver);

  // Compute the minimum and retrieve the result.
  solver_t::result_t res = solver.minimum ();

  // Display solver information.
  std::cout << solver << std::endl;

  // Process the result
  PROCESS_RESULT_UNCONSTRAINED();
}

BOOST_AUTO_TEST_SUITE_END ()
