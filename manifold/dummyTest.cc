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

#include <roboptim/core/linear-function.hh>
#include <roboptim/core/differentiable-function.hh>
#include <roboptim/core/plugin/pgsolver/pgsolver.hh>
#include <roboptim/core/plugin/pgsolver/converted-problem.hh>

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
//      (FunctionOnManifold<GenericDifferentiableFunction< T > >);
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
struct G : public GenericLinearFunction<T>
{
  //ROBOPTIM_LINEAR_FUNCTION_FWD_TYPEDEFS_(GenericLinearFunction<T>);
  ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
    (GenericLinearFunction<T>);


  G () : GenericLinearFunction<T> (22, 10, "f_n (x) = n * x")
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


BOOST_AUTO_TEST_CASE_TEMPLATE (DummyTest, T, functionTypes_t)
{
  typedef F<T> Func;

  DESC_MANIFOLD(R3, REAL_SPACE(3));
  NAMED_FUNCTION_BINDING(F_On_R3, Func, R3);

  pgs::RealSpace pos(3);pos.name() = "position";

  boost::shared_ptr<F_On_R3>
    descWrapPtr(new F_On_R3());

  Instance_F_On_R3 instWrap(descWrapPtr, pos, pos);

  solver_t::problem_t problem (instWrap);

  // Initialize solver.
  SolverFactory<solver_t> factory ("pgsolver", problem);
  solver_t& solver = factory ();

  // Solve
  solver.solve();

  std::cout << "HelloWorld" << std::endl;
  BOOST_CHECK_EQUAL(2, 2);
}


BOOST_AUTO_TEST_CASE_TEMPLATE (ConversionTest, T, functionTypes_t)
{
  typedef F<T> Func;
  typedef G<T> Gunc;

  DESC_MANIFOLD(R22, REAL_SPACE(3), roboptim::SO3, REAL_SPACE(10));
  NAMED_FUNCTION_BINDING(F_On_R22, Func, R22);
  NAMED_FUNCTION_BINDING(G_On_R22, Gunc, R22);

  pgs::RealSpace pos(3);
  pgs::SO3<pgs::ExpMapMatrix> ori;
  pgs::RealSpace joints(10);
  pgs::CartesianProduct freeFlyer(pos, ori);
  pgs::CartesianProduct robot(freeFlyer, joints);

  boost::shared_ptr<F_On_R22>
    descWrapPtr(new F_On_R22());
  boost::shared_ptr<G_On_R22>
    descWrapPtr1(new G_On_R22());

  boost::shared_ptr<Instance_F_On_R22> instWrapPtr (new Instance_F_On_R22(descWrapPtr, robot, robot));
  boost::shared_ptr<Instance_G_On_R22> instWrapPtr1 (new Instance_G_On_R22(descWrapPtr1, robot, robot));

  solver_t::problem_t problem (*instWrapPtr);

  typename Func::intervals_t bounds;
  solver_t::problem_t::scales_t scales;

  for(int i = 0; i < instWrapPtr->outputSize(); ++i) {
    bounds.push_back(Function::makeLowerInterval (25.));
    scales.push_back (1.);
  }

  problem.addConstraint
    (instWrapPtr,
     bounds, scales);

  bounds.clear();
  scales.clear();

  for(int i = 0; i < instWrapPtr->outputSize(); ++i) {
    bounds.push_back(Function::makeLowerInterval (14.));
    scales.push_back (1.);
  }

  problem.addConstraint
    (instWrapPtr1,
     bounds, scales);

  std::cout << "problem.argumentBounds().size(): " << problem.argumentBounds().size() << std::endl;

  for(size_t i = 0; i < 22; ++i)
    {
      problem.argumentBounds()[i] = std::make_pair(1l, 1l);
    }

  std::cout << "problem.argumentBounds().size(): " << problem.argumentBounds().size() << std::endl;

  roboptim::pgsolver::ConvertedProblem<solver_t>* cP = roboptim::pgsolver::ConvertedProblem<solver_t>::convertProblem(problem, robot);

  std::cout << "cP->numberOfCstr(): " << cP->numberOfCstr() << std::endl;

  Eigen::VectorXd store = Eigen::VectorXd::Zero(instWrapPtr->outputSize());
  Eigen::VectorXd store2 = Eigen::VectorXd::Zero(16);

  cP->getNonLinCstrLB(store, 0);

  std::cout << "store: " << store << std::endl;

  cP->getTangentLB(store2);

  std::cout << "store2: " << store2 << std::endl;


  std::cout << "OHAI" << std::endl;
  BOOST_CHECK_EQUAL(2, 2);
}
//*/
BOOST_AUTO_TEST_SUITE_END ()
