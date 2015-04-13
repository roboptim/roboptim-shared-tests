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

#include <roboptim/core/linear-function.hh>
#include <roboptim/core/differentiable-function.hh>
#include <roboptim/core/plugin/pgsolver/pgsolver.hh>
#include <roboptim/core/plugin/pgsolver/converted-problem.hh>

#include <manifolds/SO3.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/S2.h>
#include <manifolds/Point.h>
#include <manifolds/utils.h>

#include <roboptim/core/manifold-map/decorator/problem-on-manifold.hh>


using namespace pgs;
using namespace Eigen;
using namespace roboptim;

typedef boost::mpl::list< ::roboptim::EigenMatrixDense/*,
			  ::roboptim::EigenMatrixSparse*/> functionTypes_t;

BOOST_FIXTURE_TEST_SUITE (manifold, TestSuiteConfiguration)

template<class T>
struct SquaredNormFunc : public GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (GenericDifferentiableFunction<T>);

  SquaredNormFunc () : GenericDifferentiableFunction<T> (3, 1, "f_n (x) = n * x")
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
struct BelongsToPlane : public GenericLinearFunction<T>
{
  ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (GenericLinearFunction<T>);

  BelongsToPlane (double a, double b, double c, double d) :
    GenericLinearFunction<T> (3, 1, "f_n (x) = n * x"),
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

  DESC_MANIFOLD(R3, REAL_SPACE(3));
  NAMED_FUNCTION_BINDING(SquaredNorm_On_R3, SquaredNormFunc<T>, R3);
  NAMED_FUNCTION_BINDING(BelongsToPlane_On_R3, BelongsToPlane<T>, R3);

  pgs::RealSpace r3(3);

  boost::shared_ptr<SquaredNorm_On_R3>
    squaredNormDesc(new SquaredNorm_On_R3());

  boost::shared_ptr<BelongsToPlane_On_R3>
    belongsToPlaneDesc(new BelongsToPlane_On_R3(1, 0.5, -2, 0.6));

  boost::shared_ptr<Instance_SquaredNorm_On_R3> objFunc(new Instance_SquaredNorm_On_R3(squaredNormDesc, r3, r3));
  boost::shared_ptr<Instance_BelongsToPlane_On_R3> belongsToPlane(new Instance_BelongsToPlane_On_R3(belongsToPlaneDesc, r3, r3));
  boost::shared_ptr<Instance_SquaredNorm_On_R3> betweenTwoSpheres(new Instance_SquaredNorm_On_R3(squaredNormDesc, r3, r3));

  // Create a Roboptim problem
  solver_t::problem_t problem (*objFunc);

  typename SquaredNormFunc<T>::intervals_t bounds;
  solver_t::problem_t::scales_t scales;

  bounds.push_back(Function::makeInterval (0., 0.));
  scales.push_back (1.);

  problem.addConstraint
    (belongsToPlane,
     bounds, scales);

  bounds.clear();
  scales.clear();

  bounds.push_back(Function::makeInterval (R1*R1, R2*R2));
  scales.push_back (1.);

  problem.addConstraint
    (betweenTwoSpheres,
     bounds, scales);

#ifndef NDEBUG
  SolverFactory<solver_t> factory ("pgsolver_d", problem);
#else
  SolverFactory<solver_t> factory ("pgsolver", problem);
#endif
  solver_t& solver = factory ();

  // Solve
  solver.solve();
}

BOOST_AUTO_TEST_SUITE_END ()
