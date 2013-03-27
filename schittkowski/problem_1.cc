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

#include <iostream>

#include <boost/mpl/list.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>

#include <roboptim/core/differentiable-function.hh>
#include <roboptim/core/io.hh>
#include <roboptim/core/solver.hh>
#include <roboptim/core/solver-factory.hh>

#ifndef SOLVER_NAME
# error "please define solver name"
#endif //! PROBLEM_TYPE

#define FORWARD_TYPEDEFS()				  \
  typedef GenericDifferentiableFunction<T> parent_t;	  \
  typedef typename parent_t::result_t result_t;		  \
  typedef typename parent_t::size_type size_type;	  \
  typedef typename parent_t::argument_t argument_t;	  \
  typedef typename parent_t::gradient_t gradient_t;	  \
  typedef typename parent_t::jacobian_t jacobian_t

typedef boost::mpl::list< ::roboptim::EigenMatrixDense> functionTypes_t;

namespace roboptim
{
  namespace schittkowski
  {
    namespace problem1
    {
      struct ExpectedResult
      {
	static const double f0;
      };
      const double ExpectedResult::f0 = 909.;

      template <typename T>
      class F : public GenericDifferentiableFunction<T>
      {
      public:
	FORWARD_TYPEDEFS ();

	explicit F () throw ();
	void
	impl_compute (result_t& result, const argument_t& x) const throw ();
	void
	impl_gradient (gradient_t& grad, const argument_t& x, size_type)
	  const throw ();
      };

      template <typename T>
      F<T>::F () throw ()
	: GenericDifferentiableFunction<T>
	  (2, 1, "100 * (x[1] - x[0]^2)^2 + (1 - x[0])^2")
      {}

      template <typename T>
      void
      F<T>::impl_compute (result_t& result, const argument_t& x)
	const throw ()
      {
	result[0] = 100 * std::pow (x[1] - std::pow (x[0], 2), 2)
	  + std::pow (1 - x[0], 2);
      }

      template <typename T>
      void
      F<T>::impl_gradient (gradient_t& grad, const argument_t& x, size_type)
	const throw ()
      {
	grad[0] = -400 * x[0] * (x[1] - std::pow (x[0], 2)) - 2 * (1 - x[0]);
	grad[1] = 200 * (x[1] - std::pow (x[0], 2));
      }
    } // end of namespace problem1.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_AUTO_TEST_CASE_TEMPLATE (schittkowski_problem1, T, functionTypes_t)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem1;

  typedef Solver<GenericDifferentiableFunction<T>,
		 boost::mpl::vector<GenericLinearFunction<T>,
				    GenericDifferentiableFunction<T> > >
  solver_t;

  // Build problem.
  F<T> f;
  typename solver_t::problem_t problem (f);

  problem.argumentBounds ()[1] = F<T>::makeLowerInterval (-1.5);

  typename F<T>::argument_t x (2);
  x << -2., 1.;
  problem.startingPoint () = x;

  BOOST_CHECK_CLOSE (f (x)[0], ExpectedResult::f0, 1e-6);

  // Initialize solver.
  SolverFactory<solver_t> factory (SOLVER_NAME, problem);
  solver_t& solver = factory ();

  // Compute the minimum and retrieve the result.
  typename solver_t::result_t res = solver.minimum ();

  // Display solver information.
  std::cout << solver << std::endl;

  // Check if the minimization has succeed.
  if (res.which () != solver_t::SOLVER_VALUE)
    {
      std::cout << "A solution should have been found. Failing..."
                << std::endl
                << boost::get<SolverError> (res).what ()
                << std::endl;
      BOOST_CHECK_EQUAL (res.which (), solver_t::SOLVER_VALUE);
      return;
    }

  // Get the result.
  Result& result = boost::get<Result> (res);

  // Display the result.
  std::cout << "A solution has been found: " << std::endl
	    << result << std::endl;
}
