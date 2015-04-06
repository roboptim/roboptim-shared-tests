#ifndef COMMON_HH_U6TKTBRJ
#define COMMON_HH_U6TKTBRJ

# include <roboptim/core/differentiable-function.hh>
# include <roboptim/core/linear-function.hh>
# include <roboptim/core/decorator/manifold-map/function-on-manifold.hh>

namespace roboptim
{
  template <typename T>
  using ManifoldLinearFunction_t = FunctionOnManifold<GenericLinearFunction< T > >;
  template <typename T>
  using ManifoldDifferentiableFunction_t = FunctionOnManifold<GenericDifferentiableFunction< T > >;
}

# include "common.hh"

#endif /* end of include guard: COMMON_HH_U6TKTBRJ */
