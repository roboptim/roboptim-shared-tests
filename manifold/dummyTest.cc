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


#include "manifold/manifold_common.hh"

#include <boost/test/unit_test.hpp>

#include <roboptim/core/differentiable-function.hh>
#include <roboptim/core/plugin/pgsolver/pgsolver.hh>
#include <roboptim/core/decorator/manifold-map/manifold-map.hh>
#include <roboptim/core/decorator/manifold-map/function-on-manifold.hh>

#include <manifolds/SO3.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/S2.h>

using namespace pgs;
using namespace Eigen;
using namespace roboptim;

typedef boost::mpl::list< ::roboptim::EigenMatrixDense/*,
			  ::roboptim::EigenMatrixSparse*/> functionTypes_t;

template<class T>
struct F : public GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (GenericDifferentiableFunction<T>);

  F () : GenericDifferentiableFunction<T> (1, 3, "f_n (x) = n * x")
  {}

  void impl_compute (result_ref res, const_argument_ref ) const
  {
    res.setZero ();
  }

  void impl_gradient (gradient_ref grad, const_argument_ref,
          size_type ) const
  {
    grad.setZero ();
  }
};

BOOST_FIXTURE_TEST_SUITE (manifold, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE_TEMPLATE (DummyTest, T, functionTypes_t)
{
  typedef F<T> Func;
  
  DESC_MANIFOLD(R3, REAL_SPACE(3));
  NAMED_FUNCTION_BINDING(F_On_R3, Func, R3);

  pgs::RealSpace pos(3);
  pos.name() = "position";

  boost::shared_ptr<F_On_R3> descWrapPtr(new F_On_R3());

  Instance_F_On_R3 instWrap(descWrapPtr, pos, pos);

  solver_t::problem_t problem (instWrap);

  // Initialize solver.
#ifdef NDEBUG
  SolverFactory<solver_t> factory ("pgsolver", problem);
#else
  SolverFactory<solver_t> factory ("pgsolver_d", problem);
#endif
  solver_t& solver = factory ();

  //Solve
  solver.solve();
}

BOOST_AUTO_TEST_SUITE_END ()

