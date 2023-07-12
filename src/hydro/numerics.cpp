#include "numerics.h"

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

double Sapphire::Hydro::minmod(const std::vector<double> &values) {
  auto [min_it, max_it] = std::minmax_element(values.begin(), values.end());
  if ((*min_it) * (*max_it) < 0.0)
    return 0.0;
  else if (std::abs(*min_it) < std::abs(*max_it))
    return *min_it;
  else
    return *max_it;
}

template <int dim>
void Sapphire::Hydro::minmod(const std::vector<Tensor<1, dim>> &values,
                             const unsigned int n,
                             Tensor<1, dim> &return_value) {
  std::vector<double> component_values(n);
  for (unsigned int d = 0; d < dim; ++d) {
    for (unsigned int i = 0; i < n; ++i)
      component_values[i] = values[i][d];
    return_value[d] = minmod(component_values);
  }
}

// explicit instantiation
template void
Sapphire::Hydro::minmod<1>(const std::vector<Tensor<1, 1>> &values,
                           const unsigned int n, Tensor<1, 1> &return_value);
template void
Sapphire::Hydro::minmod<2>(const std::vector<Tensor<1, 2>> &values,
                           const unsigned int n, Tensor<1, 2> &return_value);
template void
Sapphire::Hydro::minmod<3>(const std::vector<Tensor<1, 3>> &values,
                           const unsigned int n, Tensor<1, 3> &return_value);
