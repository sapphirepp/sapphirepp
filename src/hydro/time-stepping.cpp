#include "time-stepping.h"

#include <deal.II/base/time_stepping.h>

#include <deal.II/matrix_free/operators.h> //TODO_HD: minimize includes

#include "config.h"

Sapphire::Hydro::LowStorageRungeKuttaIntegrator::LowStorageRungeKuttaIntegrator(
  const LowStorageRungeKuttaScheme scheme)
{
  TimeStepping::runge_kutta_method lsrk;
  switch (scheme)
    {
      case stage_3_order_3:
        {
          lsrk = TimeStepping::LOW_STORAGE_RK_STAGE3_ORDER3;
          break;
        }

      case stage_5_order_4:
        {
          lsrk = TimeStepping::LOW_STORAGE_RK_STAGE5_ORDER4;
          break;
        }

      case stage_7_order_4:
        {
          lsrk = TimeStepping::LOW_STORAGE_RK_STAGE7_ORDER4;
          break;
        }

      case stage_9_order_5:
        {
          lsrk = TimeStepping::LOW_STORAGE_RK_STAGE9_ORDER5;
          break;
        }

      default:
        AssertThrow(false, ExcNotImplemented());
    }
  TimeStepping::LowStorageRungeKutta<LinearAlgebra::distributed::Vector<Number>>
    rk_integrator(lsrk);
  rk_integrator.get_coefficients(ai, bi, ci);
}



unsigned int
Sapphire::Hydro::LowStorageRungeKuttaIntegrator::n_stages() const
{
  return bi.size();
}
