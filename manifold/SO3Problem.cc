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

#include <manifolds/SO3.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/S2.h>
#include <manifolds/Point.h>
#include <manifolds/utils.h>

#include <roboptim/core/manifold-map/decorator/problem-on-manifold.hh>
#include <roboptim/core/manifold-map/decorator/manifold-problem-factory.hh>


using namespace Eigen;

typedef boost::mpl::list< ::roboptim::EigenMatrixDense/*,
			  ::roboptim::EigenMatrixSparse*/> functionTypes_t;

BOOST_FIXTURE_TEST_SUITE (manifold, TestSuiteConfiguration)


size_t nPoints;
mnf::SO3<mnf::ExpMapMatrix> SO3_;
std::vector<Eigen::Vector3d> PCI; //Initial pointCloud
std::vector<Eigen::Vector3d> PCTMP; //Tmp pointCloud
std::vector<Eigen::Vector3d> PCG; //Goal pointCloud
Eigen::Matrix3d goalRot_; // Goal rotation
typedef Eigen::Map<const Eigen::Matrix3d> toMat3;

template<class T>
struct PointCloudDistFunc : public roboptim::GenericDifferentiableFunction<T>
{
  ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (roboptim::GenericDifferentiableFunction<T>);

  PointCloudDistFunc () : roboptim::GenericDifferentiableFunction<T> (9, 1, "Objective function")
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
struct RemoveOneRotation : public roboptim::GenericLinearFunction<T>
{
  ROBOPTIM_TWICE_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
  (roboptim::GenericLinearFunction<T>);

  RemoveOneRotation () :
    roboptim::GenericLinearFunction<T> (9, 2, "RemoveOneRotation")
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

  mnf::Point goalRot = SO3_.getZero();
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

  ROBOPTIM_DESC_MANIFOLD(RotSpace, roboptim::SO3);
  ROBOPTIM_NAMED_FUNCTION_BINDING(PC_Dist_On_RotSpace, PointCloudDistFunc<T>, RotSpace);
  ROBOPTIM_NAMED_FUNCTION_BINDING(Remove_Rotation_On_RotSpace, RemoveOneRotation<T>, RotSpace);
  PC_Dist_On_RotSpace pcDistDesc;
  Remove_Rotation_On_RotSpace remRotDesc;

  roboptim::ManifoldProblemFactory<solver_t::problem_t> problemFactory;

  typename RemoveOneRotation<T>::intervals_t bounds;
  bounds.push_back(roboptim::Function::makeInterval (0., 0.));
  bounds.push_back(roboptim::Function::makeInterval (0., 0.));

  problemFactory.addConstraint(remRotDesc, SO3_).setBounds(bounds);
  problemFactory.setObjective(pcDistDesc, SO3_);
  roboptim::ProblemOnManifold<solver_t::problem_t>* problem = problemFactory.getProblem();

#ifndef NDEBUG
  roboptim::SolverFactory<solver_t> factory ("pgsolver_d", *problem);
#else
  roboptim::SolverFactory<solver_t> factory ("pgsolver", *problem);
#endif
  solver_t& solver = factory ();
  // Solve
  solver.solve();

}

BOOST_AUTO_TEST_SUITE_END ()
