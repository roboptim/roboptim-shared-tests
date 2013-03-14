
// Copyright (c) 2011 CNRS
// Authors: Florent Lamiraux


// This file is part of roboptim-core-plugin-cminpack
// roboptim-core-plugin-cminpack is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.

// roboptim-core-plugin-cminpack is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// roboptim-core-plugin-cminpack  If not, see
// <http://www.gnu.org/licenses/>.

#include "shared-tests/common.hh"

namespace roboptim {
  namespace cminpack {
    using roboptim::DifferentiableFunction;
    /// Distance between a point on unit sphere and another point in R^3
    struct F : public DifferentiableFunction
    {
      F () : DifferentiableFunction
	     (2, 3,
	      "vector between unit sphere and point (x,y,z)"),
	     point_(3)
      {
	sphericalCoordinates (point_, -1.5, -1.2);
	point_*=2.;
	std::cout << "point_ = " << point_ << std::endl;
      }

      ~F () throw () {}

      void impl_compute (result_t& result, const argument_t& x) const throw ()
      {
	result.setZero ();
	double theta = x[0];
	double phi = x[1];
	sphericalCoordinates (result, theta, phi);
	result -= point_;
      }

      void impl_gradient(gradient_t& gradient, const argument_t& x,
			 size_type functionId=0) const throw ()
      {
	double theta = x[0];
	double phi = x[1];
	switch (functionId) {
	case 0:
	  gradient[0] = -sin(theta) * cos(phi);
	  gradient[1] = -cos(theta) * sin(phi);
	  break;
	case 1:
	  gradient[0] = cos(theta) * cos(phi);
	  gradient[1] = -sin(theta) * sin(phi);
	  break;
	case 2:
	  gradient[0] = 0.;
	  gradient[1] = cos(phi);
	  break;
	default:
	  abort();
	}
      }

      static void sphericalCoordinates (result_t& res, double theta, double phi)
      {
	res (0) = cos(theta) * cos(phi);
	res (1) = sin(theta) * cos(phi);
	res (2) = sin(phi);
      }
      result_t point_;
    };
    template <typename T> void initialize_problem (T& pb)
    {
      // Set initial guess (theta, phi)
      DifferentiableFunction::argument_t initialGuess(2);
      initialGuess(0) = 0.;
      initialGuess(1) = 0.;
      pb.startingPoint () = initialGuess;
    }
  } // namespace cminpack
} // namespace roboptim

