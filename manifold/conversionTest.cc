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
#include <roboptim/core/plugin/pgsolver/pgsolver.hh>
#include <roboptim/core/plugin/pgsolver/converted-problem.hh>
#include <roboptim/core/manifold-map/decorator/problem-factory.hh>
#include <roboptim/core/manifold-map/decorator/problem-on-manifold.hh>

typedef boost::mpl::list< ::roboptim::EigenMatrixDense/*,
			  ::roboptim::EigenMatrixSparse*/> functionTypes_t;

using namespace pgs;
using namespace Eigen;

template<class T>
struct F : public roboptim::GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (roboptim::GenericDifferentiableFunction<T>);

  F () : roboptim::GenericDifferentiableFunction<T> (22, 10, "f_n (x) = n * x")
  {}

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    res.setZero ();
    for (size_type i = 0; i < this->outputSize (); ++i)
      for (size_type j = 0; j < 3; ++j)
  {
    res[i] += (value_type)i * argument[0];
  }
  }

  void impl_gradient (gradient_ref grad, const_argument_ref,
          size_type functionId) const
  {
    grad.setZero ();
    for (size_type j = 0; j < 3; ++j)
      {
  grad[0] += (value_type)functionId;
      }
  }
};

template<class T>
struct G : public roboptim::GenericLinearFunction<T>
{
  ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
    (roboptim::GenericLinearFunction<T>);


  G () : roboptim::GenericLinearFunction<T> (22, 10, "f_n (x) = n * x")
  {}

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    res.setZero ();
    for (size_type i = 0; i < this->outputSize (); ++i)
      for (size_type j = 0; j < 3; ++j)
  {
    res[i] += (value_type)i * argument[0];
  }
  }

  void impl_gradient (gradient_ref grad, const_argument_ref,
          size_type functionId) const
  {
    grad.setZero ();
    for (size_type j = 0; j < 3; ++j)
      {
  grad[0] += (value_type)functionId;
      }
  }
};

BOOST_FIXTURE_TEST_SUITE (manifold, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE_TEMPLATE (ConversionTest, T, functionTypes_t)
{

  typedef ::F<T> Func;
  typedef G<T> Gunc;

  DESC_MANIFOLD(R22, REAL_SPACE(3), roboptim::SO3, REAL_SPACE(10));
  NAMED_FUNCTION_BINDING(F_On_R22, Func, R22);
  NAMED_FUNCTION_BINDING(G_On_R22, Gunc, R22);

  pgs::RealSpace pos(3);
  pgs::SO3<pgs::ExpMapMatrix> ori;
  pgs::RealSpace joints(10);
  pgs::CartesianProduct freeFlyer(pos, ori);
  pgs::CartesianProduct robot(freeFlyer, joints);

  F_On_R22 descWrapPtr;
  G_On_R22 descWrapPtr1;

  roboptim::ProblemFactory<solver_t::problem_t> probFactory;
  probFactory.setObjective(descWrapPtr, robot);
  probFactory.addConstraint(descWrapPtr1, robot);

  typename Func::intervals_t bounds;
  solver_t::problem_t::scales_t scales;

  for(int i = 0; i < descWrapPtr.fct().outputSize(); ++i) {
    bounds.push_back(roboptim::Function::makeLowerInterval (25.));
  }
  probFactory.addConstraint(descWrapPtr, robot).setBounds(bounds);

  bounds.clear();
  for(int i = 0; i < descWrapPtr.fct().outputSize(); ++i) {
    bounds.push_back(roboptim::Function::makeLowerInterval (14.));
  }
  probFactory.addConstraint(descWrapPtr1, robot).setBounds(bounds);

  solver_t::problem_t & problem = *(probFactory.getProblem());

  std::cout << "problem.argumentBounds().size(): " << problem.argumentBounds().size() << std::endl;

  for(size_t i = 0; i < 22; ++i)
    {
      problem.argumentBounds()[i] = std::make_pair(1l, 1l);
    }

  std::cout << "problem.argumentBounds().size(): " << problem.argumentBounds().size() << std::endl;


  roboptim::pgsolver::ConvertedProblem<solver_t>* cP = roboptim::pgsolver::ConvertedProblem<solver_t>::convertProblem(problem, robot);


  if (cP != nullptr)
    {
      std::cout << "cP->numberOfCstr(): " << cP->numberOfCstr() << std::endl;
    }
  Eigen::VectorXd store = Eigen::VectorXd::Zero(descWrapPtr.fct().outputSize());
  Eigen::VectorXd store2 = Eigen::VectorXd::Zero(16);

  cP->getNonLinCstrLB(store, 0);

  std::cout << "store: " << store << std::endl;

  cP->getTangentLB(store2);

  std::cout << "store2: " << store2 << std::endl;

  delete cP;

  std::cout << "OHAI" << std::endl;
  BOOST_CHECK_EQUAL(2, 2);
}

BOOST_AUTO_TEST_SUITE_END ()
