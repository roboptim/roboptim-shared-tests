// Copyright (C) 2016 by Benjamin Chrétien, CNRS-AIST JRL.
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

#include <roboptim/core/numeric-quadratic-function.hh>
#include <roboptim/core/result.hh>
#include <roboptim/core/solver-factory.hh>
#include <roboptim/core/util.hh>

namespace roboptim
{
  namespace common
  {
    namespace starting_point
    {
      template <typename T>
      struct F : public GenericNumericQuadraticFunction<T>
      {
        ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
        (GenericNumericQuadraticFunction<T>);

        explicit F () : GenericNumericQuadraticFunction<T>
                        (matrix_t (3, 3),
                         vector_t::Zero (3),
                         vector_t::Zero (1))
        {
          initialize ();
        }

        void initialize ();

        ~F ()
        {}
      };

      template <>
      void F<EigenMatrixSparse>::initialize ()
      {
        // Fill matrix A.
        Eigen::MatrixXd denseA (3, 3);
        denseA <<  2.,  0.,  0.,
                   0.,  2.,  0.,
                   0.,  0.,  2.;
        denseA *= 0.5;
        this->A () = denseA.sparseView ();
      }

      template <typename T>
      void F<T>::initialize ()
      {
        // Fill matrix A.
        this->A () <<  2.,  0.,  0.,
                       0.,  2.,  0.,
                       0.,  0.,  2.;
        this->A () *= 0.5;
      }
    } // end of namespace starting_point
  } // end of namespace common
} // end of namespace roboptim

BOOST_FIXTURE_TEST_SUITE (common, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (starting_point)
{
  using namespace roboptim;
  using namespace roboptim::common::starting_point;

  // Build problem.
  boost::shared_ptr<F<functionType_t> > f (new F<functionType_t> ());
  solver_t::problem_t problem (f);

  // Bounds on x
  for (size_t i = 0; i < 3; ++i)
    problem.argumentBounds ()[i] = Function::makeInterval (-2., 10.);

  Result resWithout (problem.function ().inputSize (),
                     problem.function ().outputSize ());
  Result resWith (problem.function ().inputSize (),
                  problem.function ().outputSize ());

  // 1) Solve without a starting point
  {
    // Initialize solver.
    SolverFactory<solver_t> factory (SOLVER_NAME, problem);
    solver_t& solver = factory ();

    // Set optimization logger
    SET_OPTIMIZATION_LOGGER (solver, "common/without-starting-point");

    // Set optional log file for debugging
    SET_LOG_FILE (solver);

    // Compute the minimum and catch the exception thrown when the minimum is
    //reached.
    solver_t::result_t res = solver.minimum ();

    // Release logger
    RELEASE_OPTIMIZATION_LOGGER ();

    // Check minimum
    BOOST_CHECK_NO_THROW (resWithout = boost::get<Result> (res));
  }

  // 2) Solve with a starting point
  {
    // Load starting point
    solver_t::problem_t::argument_t x (3);
    x << 0.5, 0.5, 0.5;
    problem.startingPoint () = x;

    // Initialize solver.
    SolverFactory<solver_t> factory (SOLVER_NAME, problem);
    solver_t& solver = factory ();

    // Set optimization logger
    SET_OPTIMIZATION_LOGGER (solver, "common/with-starting-point");

    // Set optional log file for debugging
    SET_LOG_FILE (solver);

    // Compute the minimum and catch the exception thrown when the minimum is
    //reached.
    solver_t::result_t res = solver.minimum ();

    // Release logger
    RELEASE_OPTIMIZATION_LOGGER ();

    // Check minimum
    BOOST_CHECK_NO_THROW (resWith = boost::get<Result> (res));
  }

  std::cout << resWithout << std::endl;
  std::cout << resWith << std::endl;

  const double tol = 1e-6;
  BOOST_CHECK (allclose (resWith.x, resWithout.x, tol, tol));
  BOOST_CHECK (allclose (resWith.value, resWithout.value, tol, tol));
  BOOST_CHECK (allclose (resWith.constraints, resWithout.constraints, tol, tol));
  BOOST_CHECK (allclose (resWith.lambda, resWithout.lambda, tol, tol));
}

BOOST_AUTO_TEST_SUITE_END ()
