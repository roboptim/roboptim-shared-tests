// Copyright (C) 2009 by Thomas Moulard, AIST, CNRS, INRIA.
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


#ifndef OPTIMIZATION_TESTS_HS071_BIS_HH
# define OPTIMIZATION_TESTS_HS071_BIS_HH
# include <utility>
# include <roboptim/core/twice-derivable-function.hh>

using namespace roboptim;

struct F : public TwiceDerivableFunction
{
  F () : TwiceDerivableFunction (4, 1, "a * d * (a + b + c) + d")
  {
  }

  void
  impl_compute (result_t& result, const argument_t& x) const throw ()
  {
    result.clear ();
    result (0) = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[3];
  }

  void
  impl_gradient (gradient_t& grad, const argument_t& x, size_type) const throw ()
  {
    grad.clear ();
    grad[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
    grad[1] = x[0] * x[3];
    grad[2] = x[0] * x[3] + 1;
    grad[3] = x[0] * (x[0] + x[1] + x[2]);
  }

  void
  impl_hessian (hessian_t& h, const argument_t& x, size_type) const throw ()
  {
    h.clear ();
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
};

struct G2 : public TwiceDerivableFunction
{
  G2 ()
    : TwiceDerivableFunction (4, 2, "a * b * c * d\na * a + b * b + c * c + d * d")
  {
  }

  void
  impl_compute (result_t& res, const argument_t& x) const throw ()
  {
    res.clear ();
    res (0) = x[0] * x[1] * x[2] * x[3];
    res (1) = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
  }

  void
  impl_gradient (gradient_t& grad, const argument_t& x, size_type s) const throw ()
  {
    grad.clear ();
    if (s == 0)
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

  void
  impl_hessian (hessian_t& h, const argument_t& x, size_type s) const throw ()
  {
    h.clear ();
    if (s == 0)
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
      h (0, 1) = 0.;
      h (0, 2) = 0.;
      h (0, 3) = 0.;

      h (1, 0) = 0.;
      h (1, 1) = 2.;
      h (1, 2) = 0.;
      h (1, 3) = 0.;

      h (2, 0) = 0.;
      h (2, 1) = 0.;
      h (2, 2) = 2.;
      h (2, 3) = 0.;

      h (3, 0) = 0.;
      h (3, 1) = 0.;
      h (3, 2) = 0.;
      h (3, 3) = 2.;
    }
  }
};

template <typename T, typename NLF>
void initialize_problem (T& pb)
{
  // Set bound for all variables.
  // 1. < x_i < 5. (x_i in [1.;5.])
  for (Function::size_type i = 0; i < pb.function ().inputSize (); ++i)
    pb.argumentBounds ()[i] = Function::makeInterval (1., 5.);

  // Add constraints.
  boost::shared_ptr<G2> g2 (new G2 ());

  Function::intervals_t bounds;
  bounds.push_back(Function::makeLowerInterval (25.));
  bounds.push_back(Function::makeInterval (40., 40.));

  pb.addConstraint (boost::static_pointer_cast<NLF> (g2), bounds);

  // Set the starting point.
  Function::vector_t start (pb.function ().inputSize ());
  start[0] = 1., start[1] = 5., start[2] = 5., start[3] = 1.;
  pb.startingPoint () = start;
}

#endif //! OPTIMIZATION_TESTS_HS071_BIS_HH

