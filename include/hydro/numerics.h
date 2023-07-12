/**
 * @file numerics.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement numerical functions to solve the hydrodynamics equations
 * @version 0.1
 * @date 2023-07-12
 *
 */

#ifndef HYDROSOLVER_NUMERICS_H
#define HYDROSOLVER_NUMERICS_H

#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <mpi.h>

namespace Sapphire {
namespace Hydro {
using namespace dealii;

enum class TimeSteppingScheme { ForwardEuler, ExplicitRK };
enum class FluxType { Central, Upwind, LaxFriedrich };
enum class SlopeLimiter {
  NoLimiter,
  CellAverage,
  LinearReconstruction,
  MinMod,
  MUSCL,
  GerneralizedSlopeLimiter
};

double minmod(const std::vector<double> &values);

template <int dim>
void minmod(const std::vector<Tensor<1, dim>> &values, const unsigned int n,
            Tensor<1, dim> &return_value);

} // namespace Hydro
} // namespace Sapphire
#endif