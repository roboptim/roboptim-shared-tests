// Copyright (C) 2014 by Thomas Moulard, AIST, CNRS.
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

#include <cstdlib>

namespace roboptim
{
  namespace benchmark
  {
    namespace problem17
    {
      // Same than problem_15
      template <typename T>
        class F : public GenericDifferentiableFunction<T>
      {
      public:
        ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
          (GenericDifferentiableFunction<T>);

        explicit F (size_type n);
        void
          impl_compute (result_ref result, const_argument_ref x) const;
        void
          impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
          const;

      private:
        size_type n_;
      };

      template <typename T>
        F<T>::F (size_type n)
        : GenericDifferentiableFunction<T>
          (2*n, 1, "100 (x₁ - x₀²)² + (1 - x₀)²"),
          n_ (n)
      {}

      template <typename T>
        void
        F<T>::impl_compute (result_ref result, const_argument_ref x)
        const
        {
          result[0] = 0.;
          for (size_type idx = 0; idx < 2*n_; idx += 2)
          {
            result[0] += 100 * (x[idx+1] - x[idx] * x[idx]) * (x[idx+1] - x[idx] * x[idx])
              + (1 - x[idx]) * (1 - x[idx]);
          }
          result[0] /= static_cast<value_type> (n_);
        }

      template <>
        void
        F<EigenMatrixSparse>::impl_gradient
        (gradient_ref grad, const_argument_ref x, size_type)
        const
        {
          for (size_type idx = 0; idx < 2*n_; idx += 2)
          {
            grad.coeffRef (idx) =
              400. * x[idx] * x[idx] * x[idx]
              - 400. * x[idx] * x[idx+1] + 2 * x[idx] - 2;
              grad.coeffRef (idx+1) = -200 * x[idx] * x[idx] + 200 * x[idx+1];
          }
          grad /= static_cast<value_type> (n_);
        }

      template <typename T>
        void
        F<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
        const
        {
          for (size_type idx = 0; idx < 2*n_; idx += 2)
          {
            grad[idx] =
              400. * x[idx] * x[idx] * x[idx]
              - 400. * x[idx] * x[idx+1] + 2 * x[idx] - 2;
              grad[idx+1] = -200. * x[idx] * x[idx] + 200. * x[idx+1];
          }
          grad /= static_cast<value_type> (n_);
        }

      template <typename T>
        class G : public GenericDifferentiableFunction<T>
      {
      public:
        ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
          (GenericDifferentiableFunction<T>);

        explicit G (size_type n, size_type i);
        void
          impl_compute (result_ref result, const_argument_ref x) const;
        void
          impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
          const;

      private:
        size_type idx_;
      };

      template <typename T>
        G<T>::G (size_type n, size_type i)
        : GenericDifferentiableFunction<T>
          (2*n, 1, "x₁² - x₀"),
          idx_ (i)
      {}

      template <typename T>
        void
        G<T>::impl_compute (result_ref result, const_argument_ref x)
        const
        {
          result[0] = x[idx_+1] * x[idx_+1] - x[idx_];
        }

      template <>
        void
        G<EigenMatrixSparse>::impl_gradient
        (gradient_ref grad, const_argument_ref x, size_type)
        const
        {
          grad.coeffRef (idx_) = -1.;
          grad.coeffRef (idx_+1) = 2 * x[idx_+1];
        }

      template <typename T>
        void
        G<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
        const
        {
          grad[idx_] = -1.;
          grad[idx_+1] = 2 * x[idx_+1];
        }


      template <typename T>
        class G2 : public GenericDifferentiableFunction<T>
      {
      public:
        ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
          (GenericDifferentiableFunction<T>);

        explicit G2 (size_type n, size_type i);
        void
          impl_compute (result_ref result, const_argument_ref x) const;
        void
          impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
          const;

      private:
        size_type idx_;
      };

      template <typename T>
        G2<T>::G2 (size_type n, size_type i)
        : GenericDifferentiableFunction<T>
          (2*n, 1, "x₀² - x₁"),
          idx_ (i)
      {}

      template <typename T>
        void
        G2<T>::impl_compute (result_ref result, const_argument_ref x)
        const
        {
          result[0] = x[idx_] * x[idx_] - x[idx_+1];
        }

      template <>
        void
        G2<EigenMatrixSparse>::impl_gradient
        (gradient_ref grad, const_argument_ref x, size_type)
        const
        {
          grad.coeffRef (idx_) = 2 * x[idx_];
          grad.coeffRef (idx_+1) = -1.;
        }

      template <typename T>
        void
        G2<T>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type)
        const
        {
          grad[idx_] = 2 * x[idx_];
          grad[idx_+1] = -1.;
        }
    } // end of namespace problem17.
  } // end of namespace benchmark.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (benchmark, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (benchmark_problem17)
{
  int argc = boost::unit_test::framework::master_test_suite ().argc;
  char** argv = boost::unit_test::framework::master_test_suite ().argv;

  using namespace roboptim;
  using namespace roboptim::benchmark::problem17;

  typedef F<functionType_t>::size_type size_type;
  typedef F<functionType_t>::argument_t argument_t;

  size_type n = 1;
  if (argc == 2)
  {
    n = std::atoi (argv[1]);
    if (n <= 0)
      exit (EXIT_FAILURE);
  }

  ExpectedResult expectedResult;
  expectedResult.f0 = 909.;
  expectedResult.x = ExpectedResult::argument_t::Zero (2*n);
  expectedResult.fx = 1.;

  // Tolerances for Boost checks.
  double f0_tol = 1e-4;
  double x_tol = 1e-4;
  double f_tol = 1e-4;

  // Build problem.
  boost::shared_ptr<F<functionType_t> > f (new F<functionType_t> (n));
  solver_t::problem_t problem (f);
  argument_t x (f->inputSize ());

  for (size_type i = 0; i < n; ++i)
  {
    problem.argumentBounds ()[static_cast<size_t> (2*i)]
      = F<functionType_t>::makeUpperInterval (1.);

    boost::shared_ptr<G<functionType_t> > g =
      boost::make_shared<G<functionType_t> > (n, 2*i);
    problem.addConstraint (g, G<functionType_t>::makeLowerInterval (0.));
    boost::shared_ptr<G2<functionType_t> > g2 =
      boost::make_shared<G2<functionType_t> > (n, 2*i);
    problem.addConstraint (g2, G2<functionType_t>::makeLowerInterval (0.));

    x[2*i] = -2.;
    x[2*i+1] = 1.;
  }

  problem.startingPoint () = x;

  BOOST_CHECK_SMALL_OR_CLOSE ((*f) (x)[0], expectedResult.f0, f0_tol);

  std::cout << f->inputSize () << std::endl;
  std::cout << problem.function ().inputSize () << std::endl;

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();
  // Set optimization logger
  SET_OPTIMIZATION_LOGGER (solver, "benchmark/benchmark-17");

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
