#ifndef COMMON_HH_U6TKTBRJ
#define COMMON_HH_U6TKTBRJ

# include <roboptim/core/differentiable-function.hh>
# include <roboptim/core/linear-function.hh>
# include <roboptim/core/manifold-map/decorator/manifold-map.hh>

namespace roboptim
{
  template <typename T>
  using ManifoldLinearFunction_t = GenericLinearFunction< T >;
  template <typename T>
  using ManifoldDifferentiableFunction_t = GenericDifferentiableFunction< T >;
}

# include "common.hh"

#endif /* end of include guard: COMMON_HH_U6TKTBRJ */
