/**
 * @file time-stepping.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement LowStorageRungeKuttaIntegrator class
 * @version 0.1
 * @date 2023-08-18
 *
 * @copyright Copyright (c) 2023
 *
 */


#ifndef UTILS_TIME_STEPPING_H
#define UTILS_TIME_STEPPING_H

#include <deal.II/base/time_stepping.h>

#include "parameter-flags.h"

namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;

    class LowStorageRungeKuttaIntegrator
    {
    public:
      LowStorageRungeKuttaIntegrator(const LowStorageRungeKuttaScheme scheme);

      unsigned int
      n_stages() const;

      template <typename VectorType, typename Operator> // TODO_HD: rm template
      void
      perform_time_step(const Operator &pde_operator,
                        const double    current_time,
                        const double    time_step,
                        VectorType     &solution,
                        VectorType     &vec_ri,
                        VectorType     &vec_ki) const
      {
        AssertDimension(ai.size() + 1, bi.size());

        pde_operator.perform_stage(current_time,
                                   bi[0] * time_step,
                                   ai[0] * time_step,
                                   solution,
                                   vec_ri,
                                   solution,
                                   vec_ri);

        for (unsigned int stage = 1; stage < bi.size(); ++stage)
          {
            const double c_i = ci[stage];
            pde_operator.perform_stage(current_time + c_i * time_step,
                                       bi[stage] * time_step,
                                       (stage == bi.size() - 1 ?
                                          0 :
                                          ai[stage] * time_step),
                                       vec_ri,
                                       vec_ki,
                                       solution,
                                       vec_ri);
          }
      }

    private:
      std::vector<double> bi;
      std::vector<double> ai;
      std::vector<double> ci;
    };

  } // namespace Hydro
} // namespace Sapphire


#endif