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
 * @file probe-location.h
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::VFP::ProbeLocation
 */

#ifndef VFP_PROBELOCATION_H
#define VFP_PROBELOCATION_H

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_remote_point_evaluation.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/vector.h>

#include <mpi.h>

#include <array>
#include <filesystem>
#include <vector>

#include "output-parameters.h"
#include "vfp-parameters.h"



namespace sapphirepp
{
  namespace VFP
  {
    using namespace dealii;

    /**
     * @brief PostProcessor unit to probe location in reduced phase space
     *        and reconstruct the phase space distribution.
     *
     * @tparam dim Dimension of the reduced phase space \f$ (\mathbf{x}, p) \f$
     */
    template <unsigned int dim>
    class ProbeLocation
    {
    public:
      /**
       * @brief Constructor
       *
       * The points for probe location are defined in
       * @ref VFPParameters::probe_location_points.
       *
       * @param vfp_parameters Parameters for the VFP equation
       * @param output_parameters Parameters for the output
       * @param lms_indices Map between system index \f$ i \f$
       *                    and spherical harmonic indices \f$ (l,m,s) \f$
       */
      ProbeLocation(
        const VFPParameters<dim>                       &vfp_parameters,
        const Utils::OutputParameters                  &output_parameters,
        const std::vector<std::array<unsigned int, 3>> &lms_indices);



      /**
       * @brief Prepare for reconstruction at user defined points
       *
       * @param triangulation @dealref{Triangulation}
       * @param mapping @dealref{Mapping}
       */
      void
      reinit(const Triangulation<dim> &triangulation,
             const Mapping<dim>       &mapping);



      /**
       * @brief Reconstruct the phase space distribution at all user defined
       *        points
       *
       * @param dof_handler @dealref{DoFHandler}
       * @param mapping @dealref{Mapping}
       * @param solution Current solution vector
       * @param time_step_number Current time step number
       * @param cur_time Current time
       */
      void
      reconstruct_all_points(const DoFHandler<dim>            &dof_handler,
                             const Mapping<dim>               &mapping,
                             const PETScWrappers::MPI::Vector &solution,
                             const unsigned int time_step_number = 0,
                             const double       cur_time         = 0.) const;



      /**
       * @brief Compute the real spherical harmonics
       *        \f$ Y_{lms}(\theta, \varphi) \f$
       *
       * @param theta_values Vector of theta values
       * @param phi_values Vector of phi values
       * @param lms_indices Map between system index \f$ i \f$
       *                    and spherical harmonic indices \f$ (l,m,s) \f$
       * @return Table<3, double> Real spherical harmonics `Y[i][theta][phi]`
       */
      static Table<3, double>
      compute_real_spherical_harmonics(
        const std::vector<double>                      &theta_values,
        const std::vector<double>                      &phi_values,
        const std::vector<std::array<unsigned int, 3>> &lms_indices);



    private:
      /** Output parameter */
      const Utils::OutputParameters output_parameters;

      /**
       * Map between system index \f$ i \f$
       * and spherical harmonic indices \f$ (l,m,s) \f$
       */
      const std::vector<std::array<unsigned int, 3>> lms_indices;
      /** Number of expansion coefficients */
      const unsigned int system_size;

      /** Postprocess to probe points? */
      const bool perform_probe_location;
      /** Theta values for phase space reconstruction */
      const std::vector<double> theta_values;
      /** Phi values for phase space reconstruction */
      const std::vector<double> phi_values;

      /** Real spherical harmonics `Y[i][theta][phi]` */
      const Table<3, double> real_spherical_harmonics;

      /**  Points in reduced phase space to to probe and reconstruct */
      std::vector<Point<dim>> probe_location_points;

      /** @dealref{RemotePointEvaluation,classUtilities_1_1MPI_1_1RemotePointEvaluation} */
      Utilities::MPI::RemotePointEvaluation<dim, dim> rpe_cache;



      /** @{ */
      /**
       * @brief Compute the phase space distribution
       *        \f$ f(\mathbf{x}, p, \theta, \varphi) \f$
       *        from the expansion coefficients.
       *
       * @param expansion_coefficients Expansion coefficients
       * @return std::vector<double> Phase space distribution
       *         `f[theta*n_phi + phi]`
       */
      std::vector<double>
      compute_phase_space_distribution(
        const std::vector<double> &expansion_coefficients) const;

      /**
       * @brief Output \f$ f(\theta, \varphi) \f$ in GNU splot format.
       *
       * @param f_values Phase space distribution
       * @param point_index Index of reconstructed point
       * @param time_step_number Current time step number
       * @param cur_time Current time
       */
      void
      output_gnu_splot_data(const std::vector<double> &f_values,
                            const unsigned int         point_index,
                            const unsigned int         time_step_number = 0,
                            const double               cur_time = 0.) const;



      /**
       * @brief Output \f$ f(p_x, p_y, p_z) \f$ in GNU splot format.
       *
       * @param f_values Phase space distribution
       * @param point_index Index of reconstructed point
       * @param time_step_number Current time step number
       * @param cur_time Current time
       */
      void
      output_gnu_splot_spherical_density_map(
        const std::vector<double> &f_values,
        const unsigned int         point_index,
        const unsigned int         time_step_number = 0,
        const double               cur_time         = 0.) const;
      /** @} */

      /**
       * @brief Output \f$ f_{lms} \f$ at points in location list.
       *
       * @param expansion_coefficients Values of the expansion coefficients
       * @param point_index Index of reconstructed point
       * @param time_step_number Current time step number
       * @param cur_time Current time
       */
      void
      output_f_lms(const std::vector<double> &expansion_coefficients,
                   const unsigned int         point_index,
                   const unsigned int         time_step_number = 0,
                   const double               cur_time         = 0.) const;
      /** @} */
    };
  } // namespace VFP
} // namespace sapphirepp

#endif
