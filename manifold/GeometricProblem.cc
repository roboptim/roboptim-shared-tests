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
// gnu Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with roboptim.  If not, see <http://www.gnu.org/licenses/>.


#include "manifold/manifold_common.hh"
#include "common.hh"

#include <boost/test/unit_test.hpp>

#include <manifolds/SO3.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/S2.h>
#include <manifolds/Point.h>
#include <manifolds/utils.h>

#include <roboptim/core/linear-function.hh>
#include <roboptim/core/differentiable-function.hh>
#include <roboptim/core/manifold-map/decorator/manifold-problem-factory.hh>
#include <roboptim/core/manifold-map/decorator/problem-on-manifold.hh>

typedef boost::mpl::list< ::roboptim::EigenMatrixDense/*,
			  ::roboptim::EigenMatrixSparse*/> functionTypes_t;

using namespace Eigen;

BOOST_FIXTURE_TEST_SUITE (manifold, TestSuiteConfiguration)

template<class T>
struct SquaredNormFunc : public roboptim::GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (roboptim::GenericDifferentiableFunction<T>);

  SquaredNormFunc () : roboptim::GenericDifferentiableFunction<T> (3, 1, "f_n (x) = n * x")
  {
  }

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    res.setZero ();
    for (size_type i = 0; i < 3; ++i)
      res[0] += argument[i] * argument[i];
  }

  void impl_gradient (gradient_ref grad, const_argument_ref argument,
          size_type) const
  {
    for(size_type i = 0; i < 3; ++i)
      {
	grad[i] = 2 * argument[i];
      }
  }

};
template<class T>
struct BelongsToPlane : public roboptim::GenericLinearFunction<T>
{
  ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (roboptim::GenericLinearFunction<T>);

  BelongsToPlane (double a, double b, double c, double d) :
    roboptim::GenericLinearFunction<T> (3, 1, "f_n (x) = n * x"),
    a_(a),
    b_(b),
    c_(c),
    d_(d)
  {
  }

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    res[0] = a_*argument[0] + b_*argument[1] + c_*argument[2] + d_;
  }

  void impl_gradient (gradient_ref grad, const_argument_ref,
          size_type) const
  {
    grad[0] = a_;
    grad[1] = b_;
    grad[2] = c_;
  }

private:
  double a_, b_, c_, d_;

};



BOOST_AUTO_TEST_CASE_TEMPLATE (GeometricProblemTest, T, functionTypes_t)
{
  double R1 = 0.5;
  double R2 = 1;

  roboptim::ManifoldProblemFactory<solver_t::problem_t> probFactory;

  ROBOPTIM_DESC_MANIFOLD(R3, ROBOPTIM_REAL_SPACE(3));
  ROBOPTIM_NAMED_FUNCTION_BINDING(SquaredNorm_On_R3, SquaredNormFunc<T>, R3);
  ROBOPTIM_NAMED_FUNCTION_BINDING(BelongsToPlane_On_R3, BelongsToPlane<T>, R3);

  mnf::RealSpace r3(3);

  SquaredNorm_On_R3 squaredNormDesc;
  BelongsToPlane_On_R3 belongsToPlaneDesc(1, 0.5, -2, 0.6);

  probFactory.setObjective(squaredNormDesc, r3);

  typename SquaredNormFunc<T>::intervals_t bounds;
  bounds.push_back(roboptim::Function::makeInterval (R1*R1, R2*R2));
  probFactory.addConstraint(squaredNormDesc, r3).setBounds(bounds);
  bounds.clear();

  bounds.push_back(roboptim::Function::makeInterval (0., 0.));
  probFactory.addConstraint(belongsToPlaneDesc, r3).setBounds(bounds);

  std::cout << *(probFactory.getProblem()) << std::endl;
#ifndef NDEBUG
  roboptim::SolverFactory<solver_t> factory ("pgsolver_d", *(probFactory.getProblem()));
#else
  roboptim::SolverFactory<solver_t> factory ("pgsolver", *(probFactory.getProblem()));
#endif
  solver_t& solver = factory ();

  // Solve
  solver.solve();
}

BOOST_AUTO_TEST_SUITE_END ()
