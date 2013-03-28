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

#include <log4cxx/basicconfigurator.h>

#include <roboptim/core/differentiable-function.hh>
#include <roboptim/core/io.hh>
#include <roboptim/core/solver.hh>
#include <roboptim/core/solver-factory.hh>

#ifndef SOLVER_NAME
# error "please define solver name"
#endif //! PROBLEM_TYPE

#ifndef PLUGIN_PATH
# error "please define plug-in path"
#endif //! PROBLEM_TYPE

#ifndef FUNCTION_TYPE
# error "please define function type"
#endif //! PROBLEM_TYPE

typedef FUNCTION_TYPE functionType_t;

#define FORWARD_TYPEDEFS()					  \
  typedef GenericTwiceDifferentiableFunction<T> parent_t;	  \
  typedef typename parent_t::result_t result_t;			  \
  typedef typename parent_t::size_type size_type;		  \
  typedef typename parent_t::argument_t argument_t;		  \
  typedef typename parent_t::gradient_t gradient_t;		  \
  typedef typename parent_t::jacobian_t jacobian_t;		  \
  typedef typename parent_t::hessian_t hessian_t

typedef boost::mpl::list< ::roboptim::EigenMatrixDense> functionTypes_t;

struct TestSuiteConfiguration
{
  TestSuiteConfiguration ()
  {
    log4cxx::BasicConfigurator::configure ();

    lt_dlinit();
    BOOST_REQUIRE_EQUAL (lt_dlsetsearchpath (PLUGIN_PATH), 0);
  }

  ~TestSuiteConfiguration ()
  {
    lt_dlexit ();
  }
};

namespace roboptim
{
  namespace schittkowski
  {
    namespace problem71b
    {
      struct ExpectedResult
      {
	static const double f0;
	static const double x[];
	static const double fx;
      };
      const double ExpectedResult::f0 = 16.;
      const double ExpectedResult::x[] = {1., 4.742994, 3.8211503, 1.3794082};
      const double ExpectedResult::fx = 17.0140173;


      template <typename T>
      struct F : public GenericTwiceDifferentiableFunction<T>
      {
	FORWARD_TYPEDEFS ();

	F ()
	  : GenericTwiceDifferentiableFunction<T>
	    (4, 1, "x1 * x4 * (x1 + x2 + x3) + x3")
	{}

	void
	impl_compute (result_t& result, const argument_t& x) const throw ()
	{
	  result[0] = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
	}

	void
	impl_gradient (gradient_t& grad, const argument_t& x, size_type)
	  const throw ();

	void
	impl_hessian (hessian_t& h, const argument_t& x, size_type)
	  const throw ();
      };

      template <typename T>
      struct G : public GenericTwiceDifferentiableFunction<T>
      {
	FORWARD_TYPEDEFS ();

	G ()
	  : GenericTwiceDifferentiableFunction<T>
	    (4, 2, "a * b * c * d\na * a + b * b + c * c + d * d")
	{
	}

	void
	impl_compute (result_t& res, const argument_t& x) const throw ()
	{
	  res.setZero ();
	  res (0) = x[0] * x[1] * x[2] * x[3];
	  res (1) = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
	}

	void
	impl_gradient (gradient_t& grad, const argument_t& x, size_type)
	  const throw ();

	void
	impl_hessian (hessian_t& h, const argument_t& x, size_type)
	  const throw ();
      };

      template <>
      void
      F<EigenMatrixSparse>::impl_gradient
      (gradient_t& grad, const argument_t& x, size_type)
	const throw ()
      {
	grad.setZero ();
	grad.insert (0) = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
	grad.insert (1) = x[0] * x[3];
	grad.insert (2) = x[0] * x[3] + 1;
	grad.insert (3) = x[0] * (x[0] + x[1] + x[2]);
      }

      template <typename T>
      void
      F<T>::impl_gradient
      (gradient_t& grad, const argument_t& x, size_type)
	const throw ()
      {
	grad.setZero ();
	grad[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
	grad[1] = x[0] * x[3];
	grad[2] = x[0] * x[3] + 1;
	grad[3] = x[0] * (x[0] + x[1] + x[2]);
      }


      template <>
      void
      G<EigenMatrixSparse>::impl_gradient
      (gradient_t& grad, const argument_t& x, size_type functionId)
	const throw ()
      {
	grad.setZero ();
	if (functionId == 0)
	  {
	    grad.insert (0) = 2 * x[0];
	    grad.insert (1) = 2 * x[1];
	    grad.insert (2) = 2 * x[2];
	    grad.insert (3) = 2 * x[3];
	  }
	else
	  {
	    grad.insert (0) = x[1] * x[2] * x[3];
	    grad.insert (1) = x[0] * x[2] * x[3];
	    grad.insert (2) = x[0] * x[1] * x[3];
	    grad.insert (3) = x[0] * x[1] * x[2];
	  }
      }

      template <typename T>
      void
      G<T>::impl_gradient
      (gradient_t& grad, const argument_t& x, size_type functionId)
	const throw ()
      {
	grad.setZero ();
	if (functionId == 0)
	  {
	    grad[0] = x[1] * x[2] * x[3];
	    grad[1] = x[0] * x[2] * x[3];
	    grad[2] = x[0] * x[1] * x[3];
	    grad[3] = x[0] * x[1] * x[2];
	  }
	else
	  {
	    grad[0] = 2 * x[0];
	    grad[1] = 2 * x[1];
	    grad[2] = 2 * x[2];
	    grad[3] = 2 * x[3];
	  }
      }

      template <>
      void
      F<EigenMatrixSparse>::impl_hessian
      (hessian_t& h, const argument_t& x, size_type) const throw ()
      {
	h.setZero ();
	h.insert (0, 0) = 2 * x[3];
	h.insert (0, 1) = x[3];
	h.insert (0, 2) = x[3];
	h.insert (0, 3) = 2 * x[0] + x[1] + x[2];

	h.insert (1, 0) = x[3];
	h.insert (1, 3) = x[0];

	h.insert (2, 0) = x[3];
	h.insert (2, 3) = x[1];

	h.insert (3, 0) = 2 * x[0] + x[1] + x[2];
	h.insert (3, 1) = x[0];
	h.insert (3, 2) = x[0];
      }

      template <typename T>
      void
      F<T>::impl_hessian (hessian_t& h, const argument_t& x, size_type)
	const throw ()
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
      G<EigenMatrixSparse>::impl_hessian
      (hessian_t& h, const argument_t& x, size_type functionId) const throw ()
      {
	if (functionId == 0)
	  {
	    h.insert (0, 1) = x[2] * x[3];
	    h.insert (0, 2) = x[1] * x[3];
	    h.insert (0, 3) = x[1] * x[2];

	    h.insert (1, 0) = x[2] * x[3];
	    h.insert (1, 2) = x[0] * x[3];
	    h.insert (1, 3) = x[0] * x[2];

	    h.insert (2, 0) = x[1] * x[3];
	    h.insert (2, 1) = x[0] * x[3];
	    h.insert (2, 3) = x[0] * x[1];

	    h.insert (3, 0) = x[1] * x[2];
	    h.insert (3, 1) = x[0] * x[2];
	    h.insert (3, 2) = x[0] * x[1];
	  }
	else
	  {
	    h.insert (0, 0) = 2.;
	    h.insert (1, 1) = 2.;
	    h.insert (2, 2) = 2.;
	    h.insert (3, 3) = 2.;
	  }
      }


      template <typename T>
      void
      G<T>::impl_hessian
      (hessian_t& h, const argument_t& x, size_type functionId)
	const throw ()
      {
	if (functionId == 0)
	  {
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
	else
	  {
	    h (0, 0) = 2.;
	    h (1, 1) = 2.;
	    h (2, 2) = 2.;
	    h (3, 3) = 2.;
	  }
      }
    } // end of namespace problem71b.
  } // end of namespace schittkowski.
} // end of namespace roboptim.

BOOST_FIXTURE_TEST_SUITE (schittkowski, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (problem_71b)
{
  using namespace roboptim;
  using namespace roboptim::schittkowski::problem71b;

  typedef Solver<
    GenericDifferentiableFunction<functionType_t>,
    boost::mpl::vector<GenericLinearFunction<functionType_t>,
		       GenericDifferentiableFunction<functionType_t> > >
    solver_t;

  // Build problem.
  F<functionType_t> f;
  typename solver_t::problem_t problem (f);

  // Set bound for all variables.
  // 1. < x_i < 5. (x_i in [1.;5.])
  for (std::size_t i = 0;
       i < static_cast<std::size_t> (problem.function ().inputSize ()); ++i)
    problem.argumentBounds ()[i] = Function::makeInterval (1., 5.);

  // Add constraints.
  boost::shared_ptr<G<functionType_t> > g (new G<functionType_t> ());

  typename F<functionType_t>::intervals_t bounds;
  bounds.push_back(Function::makeLowerInterval (25.));
  bounds.push_back(Function::makeInterval (40., 40.));

  typename solver_t::problem_t::scales_t scales;
  scales.push_back (1.);
  scales.push_back (1.);

  problem.addConstraint
    (boost::static_pointer_cast<
      GenericDifferentiableFunction<functionType_t>  > (g),
     bounds, scales);

  // Set the starting point.
  typename F<functionType_t>::argument_t x (4);
  x << 1., 5., 5., 1.;
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

  // Check final x.
  for (unsigned i = 0; i < result.x.size (); ++i)
    BOOST_CHECK_CLOSE (result.x[i], ExpectedResult::x[i], 1e-3);

  // Check final value.
  BOOST_CHECK_CLOSE (1. + result.value[0], 1. + ExpectedResult::fx, 1e-3);

  // Display the result.
  std::cout << "A solution has been found: " << std::endl
  	    << result << std::endl;
}

BOOST_AUTO_TEST_SUITE_END ()
