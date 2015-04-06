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

#include <manifolds/SO3.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/S2.h>

using namespace pgs;
using namespace Eigen;
using namespace roboptim;

//namespace roboptim
//{
//  namespace manifold
//  {
//    template <typename T>
//    struct F : public GenericDifferentiableFunction< T >
//    {
//      ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
//      (InstanceWrapper<GenericDifferentiableFunction< T > >);
//      F () : GenericDifferentiableFunction<T> (3, 1, "f_n (x) = empty")
//      {}
//    };
//  } // end of namespace manifold
//} // end of namespace roboptim

typedef boost::mpl::list< ::roboptim::EigenMatrixDense/*,
			  ::roboptim::EigenMatrixSparse*/> functionTypes_t;

template<class T>
struct F : public GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (GenericDifferentiableFunction<T>);

  F () : GenericDifferentiableFunction<T> (22, 10, "f_n (x) = n * x")
  {}

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    res.setZero ();
    for (size_type i = 0; i < this->outputSize (); ++i)
      for (size_type j = 0; j < 3; ++j)
  {
    res[i] += (value_type)i * argument[19 + j];
  }
  }

  void impl_gradient (gradient_ref grad, const_argument_ref,
          size_type functionId) const
  {
    grad.setZero ();
    for (size_type j = 0; j < 3; ++j)
      {
  grad[19 + j] += (value_type)functionId;
      }
  }
};

BOOST_FIXTURE_TEST_SUITE (manifold, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE_TEMPLATE (DummyTest, T, functionTypes_t)
{
  typedef F<T> Func;

  DESC_MANIFOLD(R3, REAL_SPACE(3));
  //NAMED_FUNCTION_BINDING(F_On_R3, Func, R3);
  roboptim::FunctionOnManifold<DifferentiableFunction> a;

  //DESC_MANIFOLD(FreeFlyerPlus10, REAL_SPACE(10), roboptim::SO3, REAL_SPACE(3));
  //NAMED_FUNCTION_BINDING(F_On_FreeFlyerPlus10, Func, FreeFlyerPlus10);

  //pgs::RealSpace pos(3);pos.name() = "position";

  //boost::shared_ptr<F_On_R3>
    //descWrapPtr(new F_On_R3());

  //Instance_F_On_R3 instWrap(descWrapPtr, pos, R3);

  //solver_t::problem_t problem (instWrap);

  //// Initialize solver.
  //SolverFactory<solver_t> factory ("pgsolver", problem);
  //solver_t& solver = factory ();

  //std::cout << solver << std::endl; 

  //std::cout << "HelloWorld" << std::endl;
  //BOOST_CHECK_EQUAL(2, 2);
}

BOOST_AUTO_TEST_SUITE_END ()

