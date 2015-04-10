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

/*
BOOST_AUTO_TEST_CASE_TEMPLATE (DummyTest, T, functionTypes_t)
{
  typedef F<T> Func;

  std::cout << "ZA DUMMY_TEST" << std::endl;

  DESC_MANIFOLD(R3, REAL_SPACE(3));
  NAMED_FUNCTION_BINDING(F_On_R3, Func, R3);

  pgs::RealSpace pos(3);pos.name() = "position";

  boost::shared_ptr<F_On_R3>
    descWrapPtr(new F_On_R3());

  Instance_F_On_R3 instWrap(descWrapPtr, pos, pos);

  solver_t::problem_t problem (instWrap);

  // Initialize solver.
  SolverFactory<solver_t> factory ("pgsolver_d", problem);
  solver_t& solver = factory ();

  // Solve
  solver.solve();

  std::cout << "HelloWorld" << std::endl;
  BOOST_CHECK_EQUAL(2, 2);
}
*/

/*
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


  if (cP != nullptr)
    {
      std::cout << "cP->numberOfCstr(): " << cP->numberOfCstr() << std::endl;
    }
  Eigen::VectorXd store = Eigen::VectorXd::Zero(instWrapPtr->outputSize());
  Eigen::VectorXd store2 = Eigen::VectorXd::Zero(16);

  cP->getNonLinCstrLB(store, 0);

  std::cout << "store: " << store << std::endl;

  cP->getTangentLB(store2);

  std::cout << "store2: " << store2 << std::endl;

  delete cP;

  std::cout << "OHAI" << std::endl;
  BOOST_CHECK_EQUAL(2, 2);
}
*/

// ---- //

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

/*
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

  SolverFactory<solver_t> factory ("pgsolver_d", problem);
  solver_t& solver = factory ();

  // Solve
  solver.solve();
}
*/

// ---- //

size_t nPoints;
pgs::SO3<pgs::ExpMapMatrix> SO3_;
std::vector<Eigen::Vector3d> PCI; //Initial pointCloud
std::vector<Eigen::Vector3d> PCTMP; //Tmp pointCloud
std::vector<Eigen::Vector3d> PCG; //Goal pointCloud
Eigen::Matrix3d goalRot_; // Goal rotation
typedef Eigen::Map<const Eigen::Matrix3d> toMat3;

template<class T>
struct PointCloudDistFunc : public GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (GenericDifferentiableFunction<T>);

  PointCloudDistFunc () : GenericDifferentiableFunction<T> (9, 1, "Objective function")
  {
    rot = Eigen::Matrix3d::Zero();
    dist = Eigen::Vector3d::Zero();
    outRepPoint = Eigen::VectorXd::Zero(9);
  }

  mutable Eigen::Matrix3d rot;
  mutable Eigen::Vector3d dist;
  mutable Eigen::VectorXd outRepPoint;

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    res.setZero ();
    dist.setZero();

    rot = toMat3(argument.data());
    double out = 0;
    for (size_t i = 0; i<nPoints; ++i)
    {
      dist = PCG[i]-rot*PCI[i];
      out += dist.transpose()*dist;
    }
    out = out / static_cast<double>(nPoints);

    res[0] = out;
  }

  void impl_gradient (gradient_ref grad, const_argument_ref,
          size_type) const
  {
    grad.setZero();
    outRepPoint.setZero();
    for( size_t i = 0; i < nPoints; ++i)
    {
      outRepPoint(0) = -2*PCG[i](0)*PCI[i](0),
      outRepPoint(1) = -2*PCG[i](1)*PCI[i](0),
      outRepPoint(2) = -2*PCG[i](2)*PCI[i](0),
      outRepPoint(3) = -2*PCG[i](0)*PCI[i](1),
      outRepPoint(4) = -2*PCG[i](1)*PCI[i](1),
      outRepPoint(5) = -2*PCG[i](2)*PCI[i](1),
      outRepPoint(6) = -2*PCG[i](0)*PCI[i](2),
      outRepPoint(7) = -2*PCG[i](1)*PCI[i](2),
      outRepPoint(8) = -2*PCG[i](2)*PCI[i](2);
      grad += outRepPoint;
    }
    grad.array() /= static_cast<double>(nPoints);
  }

};

template<class T>
struct RemoveOneRotation : public GenericLinearFunction<T>
{
  ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (GenericLinearFunction<T>);

  RemoveOneRotation () :
    GenericLinearFunction<T> (9, 2, "RemoveOneRotation")
  {
    rot = Eigen::Matrix3d::Zero();
  }

  mutable Eigen::Matrix3d rot;

  void impl_compute (result_ref res, const_argument_ref argument) const
  {
    rot = toMat3(argument.data());
    res[0] = rot(2, 0);
    res[1] = rot(2, 1);
  }

  void impl_gradient (gradient_ref grad, const_argument_ref,
          size_type functionId) const
  {
    if (functionId == 0)
      {
	grad << 0, 0, 1, 0, 0, 0, 0, 0, 0;
      }
    else
      {
	grad << 0, 0, 0, 0, 0, 1, 0, 0, 0;
      }
  }

};

BOOST_AUTO_TEST_CASE_TEMPLATE (SO3ProblemTest, T, functionTypes_t)
{
  Eigen::Vector3d v;
  v << 0.10477, 0.03291, -0.19174;

  pgs::Point goalRot = SO3_.getZero();
  goalRot.increment(v);
  std::cout << "goalRot = \n" << goalRot << std::endl;
  goalRot_ = toMat3(goalRot[0].data());
  nPoints = 300;
  PCI.resize(nPoints);
  PCG.resize(nPoints);
  double nPoints_d = static_cast<double>(nPoints);
  for (size_t i = 0; i < nPoints/3; ++i)
    PCI[i] << cos(3*static_cast<double>(i)*2*M_PI/nPoints_d), sin(3*static_cast<double>(i)*2*M_PI/nPoints_d), 0;
  for (size_t i = nPoints/3; i < 2*nPoints/3; ++i)
    PCI[i] << 0, cos(3*static_cast<double>(i)*2*M_PI/nPoints_d), sin(3*static_cast<double>(i)*2*M_PI/nPoints_d);
  for (size_t i = 2*nPoints/3; i < nPoints; ++i)
    PCI[i] << cos(3*static_cast<double>(i)*2*M_PI/nPoints_d), 0, sin(3*static_cast<double>(i)*2*M_PI/nPoints_d);

  for (size_t i = 0; i<nPoints; ++i)
    {
      PCG[i] = goalRot_*PCI[i];
    }

  DESC_MANIFOLD(RotSpace, roboptim::SO3);
  NAMED_FUNCTION_BINDING(PC_Dist_On_RotSpace, PointCloudDistFunc<T>, RotSpace);
  NAMED_FUNCTION_BINDING(Remove_Rotation_On_RotSpace, RemoveOneRotation<T>, RotSpace);

  boost::shared_ptr<PC_Dist_On_RotSpace>
    pcDistDesc(new PC_Dist_On_RotSpace());

  boost::shared_ptr<Remove_Rotation_On_RotSpace>
    remRotDesc(new Remove_Rotation_On_RotSpace());

  Instance_PC_Dist_On_RotSpace objFunc(pcDistDesc, SO3_, SO3_);

  boost::shared_ptr<Instance_Remove_Rotation_On_RotSpace> remRotCnstr(new Instance_Remove_Rotation_On_RotSpace(remRotDesc, SO3_, SO3_));

  roboptim::ProblemOnManifold<solver_t::problem_t> problem(SO3_, objFunc);

  typename RemoveOneRotation<T>::intervals_t bounds;
  solver_t::problem_t::scales_t scales;

  bounds.push_back(Function::makeInterval (0., 0.));
  bounds.push_back(Function::makeInterval (0., 0.));
  scales.push_back (1.);
  scales.push_back (1.);

  problem.addConstraint
    (remRotCnstr,
     bounds, scales);


  SolverFactory<solver_t> factory ("pgsolver_d", problem);
  solver_t& solver = factory ();

  // Solve
  solver.solve();

}

BOOST_AUTO_TEST_SUITE_END ()
