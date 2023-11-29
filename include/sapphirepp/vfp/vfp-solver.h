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

/**
 * @file vfp-solver.h
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::VFP::VFPSolver
 */

#ifndef VFP_VFPSOLVER_H
#define VFP_VFPSOLVER_H

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/types.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/vector_tools.h>

#include <mpi.h>

#include <type_traits>

#include "config.h"
#include "output-parameters.h"
#include "pde-system.h"
#include "upwind-flux.h"
#include "vfp-flags.h"
#include "vfp-parameters.h"



namespace sapphirepp
{
  namespace VFP
  {
    using namespace dealii;


    /**
     * @brief This class solves the Vlasov-Fokker-Planck equation
     *
     * Solve the Vlasov-Fokker-Planck equation specified with by the @ref
     * VFPFlags. We distinguish between two different dimensions of the problem:
     *
     *  - The dimension of the \f$ (\mathbf{x}) \f$-space, called configuration
     *    space (`dim_cs`). Depending on the setup, this can be 1, 2 or 3.
     *  - The total dimension of the problem in the reduced phase space \f$
     *    (\mathbf{x}, p) \f$, `dim_ps`. In the momentum depended case this
     *    equals `dim_ps = dim_cs + 1`, in the momentum independent
     *    (transport-only) case `dim_ps = dim_cs`.
     *
     * Since we decompose the \f$ \mathbf{p} \f$ momentum part into spherical
     * harmonics, we always solve for the full 3 dimensions in momentum space.
     *
     *
     * @note Due to limitations of @dealii, the total dimension of the problem
     *       must be smaller or equal to 3, `dim_ps <= 3`.
     *
     *
     * @tparam dim Total dimension of the problem in reduced phase space \f$
     *         (\mathbf{x}, p) \f$ (`dim_ps`)
     */
    template <unsigned int dim>
    class VFPSolver
    {
    private:
      /** Is the momentum term activated? */
      static constexpr bool momentum =
        (vfp_flags & VFPFlags::momentum) != VFPFlags::none ? true : false;
      /** Dimension in reduced phase space */
      static constexpr int dim_ps = dim;
      /** Dimension of the configuration space */
      static constexpr int dim_cs = dim - momentum;
      /** Do we use a logarithmic momentum variable? */
      static constexpr bool logarithmic_p =
        (vfp_flags & VFPFlags::linear_p) == VFPFlags::none;

      /**
       * Type of triangulation.
       *
       * @note `parallel::distributed:Triangulation` does not allow 1D. To
       *       maintain the possibility to compute 1D scenarios, we replace it
       *       with `parallel::shared::Triangulation`.
       */
      using Triangulation = typename std::conditional<
        dim_ps != 1,
        parallel::distributed::Triangulation<dim_ps>,
        parallel::shared::Triangulation<dim_ps>>::type;



    public:
      VFPSolver(const VFPParameters<dim_ps>   &vfp_parameters,
                const PhysicalParameters      &physical_parameters,
                const Utils::OutputParameters &output_parameters);



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

      MPI_Comm           mpi_communicator;
      const unsigned int n_mpi_procs;
      const unsigned int rank;

      ConditionalOStream      pcout;
      Utils::OutputParameters output_parameters;

      Triangulation triangulation;

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

      const unsigned int expansion_order;
      const unsigned int num_exp_coefficients;

      TimerOutput timer;



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
    };

  } // namespace VFP
} // namespace sapphirepp
#endif
