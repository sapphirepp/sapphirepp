// -----------------------------------------------------------------------------
//
// Copyright (C) 2023 by the Sapphire++ authors
//
// This file is part of Sapphire++.
//
// Sapphire++ is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// Sapphire++ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Sapphire++. If not, see <https://www.gnu.org/licenses/>.
//
// -----------------------------------------------------------------------------

#ifndef VFP_VFPSOLVER_H
#define VFP_VFPSOLVER_H

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/fe/fe_update_flags.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h> // neeed for periodic boundary conditions
                                     // (collect_periodic_faces) and
                                     // partitioning of the distributed mesh
// #include <deal.II/grid/grid_refinement.h>
// #include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
// #include <deal.II/distributed/grid_refinement.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
// PetscWrappers
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparsity_tools.h> // need for distribute_sparsity_pattern
                                        // to fit into a global numbering
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/lac/vector_operation.h>

#include <deal.II/meshworker/assemble_flags.h>
#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_project.h>

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "config.h"
#include "output-parameters.h"
#include "particle-functions.h"
#include "pde-system.h"
#include "sapphirepp-logstream.h"
#include "upwind-flux.h"
#include "vfp-flags.h"
#include "vfp-parameters.h"



namespace Sapphire
{
  namespace VFP
  {
    using namespace dealii;

    template <int dim>
    class VFPSolver
    {
    private:
      static constexpr bool logarithmic_p =
        (vfp_flags & VFPFlags::linear_p) == VFPFlags::none;
      static constexpr bool momentum =
        (vfp_flags & VFPFlags::momentum) != VFPFlags::none ? true : false;
      static constexpr int dim_ps = dim;
      static constexpr int dim_cs = dim - momentum;

    public:
      VFPSolver(const VFPParameters<dim_ps>           &vfp_parameters,
                const PhysicalParameters              &physical_parameters,
                const Utils::OutputParameters<dim_ps> &output_parameters);
      void
      run();

      double
      compute_global_error(
        const Function<dim_ps>         &exact_solution,
        const VectorTools::NormType    &cell_norm,
        const VectorTools::NormType    &global_norm,
        const Function<dim_ps, double> *weight = nullptr) const;

      unsigned int
      get_n_dofs() const;

    private:
      const VFPParameters<dim_ps> vfp_parameters;
      const PhysicalParameters    physical_parameters;

      // Triangulation
      void
      make_grid();
      // Setup data structures for the linear system
      void
      setup_system();
      // Matrix assembly
      void
      assemble_mass_matrix();
      void
      assemble_dg_matrix(const double time);
      // Time stepping methods
      void
      theta_method(const double time, const double time_step);
      void
      explicit_runge_kutta(const double time, const double time_step);
      void
      low_storage_explicit_runge_kutta(const double time,
                                       const double time_step);
      // Output
      void
      output_results(const unsigned int time_step_number);

      // auxiliary functions
      void
      project(const Function<dim>        &f,
              PETScWrappers::MPI::Vector &projected_function);
      // compute the source term
      void
      compute_source_term(const Function<dim> &source_function);

      MPI_Comm           mpi_communicator;
      const unsigned int n_mpi_procs;
      const unsigned int rank;

      ConditionalOStream              pcout;
      Utils::OutputParameters<dim_ps> output_parameters;

      // NOTE: parallel::distributed:Triangulation does not allow 1D. This
      // excludes the 1D transport (i.e. no momentum terms) only case. But I
      // would like to maintain the possibility to compute 1D scenarios. In
      // Step-17, Step-18 it is explained how to deal with copys of the
      // triangulation and the dof handler on every mpi process: It is enough to
      // replace parallel::distributed:Triangulation with
      // parallel::shared::Triangulation. We use a bit of C++ metaprogramming
      // tricky to decide which triangulation to use.
      typename std::conditional<dim_ps != 1,
                                parallel::distributed::Triangulation<dim_ps>,
                                parallel::shared::Triangulation<dim_ps>>::type
        triangulation;

      DoFHandler<dim_ps> dof_handler;

      IndexSet locally_owned_dofs;
      IndexSet locally_relevant_dofs;

      // NOTE: The explicit use of a mapping is most likely related to the usage
      // of mesh_loop as well
      const MappingQ1<dim_ps> mapping;
      const FESystem<dim_ps>  fe; // TODO: const is probably wrong

      // NOTE: Quadratures are members of this class (and not e.g. part of the
      // assemble_system method), because I am using the mesh_loop function
      const QGauss<dim_ps>     quadrature;
      const QGauss<dim_ps - 1> quadrature_face;

      // The constraints object is used for the distribute_local_to_global()
      // function
      const AffineConstraints<double> constraints;

      // PDE System
      PDESystem                                   pde_system;
      UpwindFlux<dim_ps, momentum, logarithmic_p> upwind_flux;

      SparsityPattern                  sparsity_pattern;
      PETScWrappers::MPI::SparseMatrix mass_matrix;
      PETScWrappers::MPI::SparseMatrix dg_matrix;
      PETScWrappers::MPI::SparseMatrix system_matrix;

      PETScWrappers::MPI::Vector system_rhs;
      PETScWrappers::MPI::Vector locally_owned_previous_solution;
      PETScWrappers::MPI::Vector locally_relevant_current_solution;

      PETScWrappers::MPI::Vector locally_owned_current_source;

      const int          expansion_order      = vfp_parameters.expansion_order;
      const unsigned int num_exp_coefficients = static_cast<unsigned int>(
        (expansion_order + 1) * (expansion_order + 1));

      TimerOutput timer;
    };

  } // namespace VFP
} // namespace Sapphire
#endif
