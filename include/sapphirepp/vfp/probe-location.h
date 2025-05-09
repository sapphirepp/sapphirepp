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
       * @brief Probe all user defined point in the reduced phase space
       *        and perform phase space reconstruction.
       *
       * @param dof_handler @dealref{DoFHandler}
       * @param mapping @dealref{Mapping}
       * @param solution Current solution vector
       * @param time_step_number Current time step number
       * @param cur_time Current time
       * @param function_symbol The symbol used for the distribution function.
       * Default, is \f$ f \f$ for an unscaled distribution function
       * and \f$ g \f$ for \f$p^3 f \f$.
       */
      void
      probe_all_points(const DoFHandler<dim>            &dof_handler,
                       const Mapping<dim>               &mapping,
                       const PETScWrappers::MPI::Vector &solution,
                       const unsigned int                time_step_number = 0,
                       const double                      cur_time         = 0.,
                       const std::string &function_symbol = "f") const;



      /**
       * @brief Compute the real spherical harmonics
       *        \f$ Y_{lms}(\theta, \varphi) \f$
       *
       * @param cos_theta_values Vector of \f$ \cos(\theta) \f$ values
       * @param phi_values Vector of \f$ \phi \f$ values
       * @param lms_indices Map between system index \f$ i \f$
       *                    and spherical harmonic indices \f$ (l,m,s) \f$
       * @return Table<3, double> Real spherical harmonics `Y[i][theta][phi]`
       */
      static Table<3, double>
      compute_real_spherical_harmonics(
        const std::vector<double>                      &cos_theta_values,
        const std::vector<double>                      &phi_values,
        const std::vector<std::array<unsigned int, 3>> &lms_indices);

      /**
       * @brief Test the phase space reconstruction.
       *
       * Call this function, for example, in the constructor of the
       * VFPSolver. It will output all spherical harmonics up to \f$
       * l_{\mathrm{max}} \f$. Results can be compared with Appendix A in
       * @cite Schween2025 .
       *
       * @note The execution of Sapphire++ will be aborted after the test is
       * executed. Moreover, the test only runs successfully, if the parameter
       * files contains the section "Probe Location" with at least one point
       * specified.
       */
      void
      test_phase_space_reconstruction();

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
      /** Perform phase space reconstruction? */
      const bool perform_phase_space_reconstruction;
      /** Cos theta values for phase space reconstruction */
      const std::vector<double> cos_theta_values;
      /** Phi values for phase space reconstruction */
      const std::vector<double> phi_values;

      /** Real spherical harmonics `Y[i][theta][phi]` */
      const Table<3, double> real_spherical_harmonics;

      /**  Points in reduced phase space to probe and reconstruct */
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
       * @brief Output \f$ f(\cos\theta, \varphi) \f$ or $\f$ g = p^3
       * f(\cos\theta, \varphi) \f$ in GNU splot format.
       *
       * @param distribution_func_values Phase space distribution
       * @param point_index Index of reconstructed point
       * @param time_step_number Current time step number
       * @param cur_time Current time
       * @param function_symbol The symbol used for the distribution function.
       * Default, is \f$ f \f$ for an unscaled distribution function
       * and \f$ g \f$ for \f$p^3 f \f$.
       */
      void
      output_gnu_splot_data(const std::vector<double> &distribution_func_values,
                            const unsigned int         point_index,
                            const unsigned int         time_step_number = 0,
                            const double               cur_time         = 0.,
                            const std::string &function_symbol = "f") const;



      /**
       * @brief Output \f$ f(p_x, p_y, p_z) \f$ or \f$ g = p^3 f(p_x, p_y, p_z)
       * \f$ in GNU splot format.
       *
       * @param distribution_func_values Phase space distribution
       * @param point_index Index of reconstructed point
       * @param time_step_number Current time step number
       * @param cur_time Current time
       * @param function_symbol The symbol used for the distribution function.
       * Default, is \f$ f \f$ for an unscaled distribution function
       * and \f$ g \f$ for \f$p^3 f \f$.
       */
      void
      output_gnu_splot_spherical_density_map(
        const std::vector<double> &distribution_func_values,
        const unsigned int         point_index,
        const unsigned int         time_step_number = 0,
        const double               cur_time         = 0.,
        const std::string         &function_symbol  = "f") const;
      /** @} */



      /**
       * @brief Outputs the values of the expansion coefficients at points in
       * location list.
       *
       * @param expansion_coefficients Values of the expansion coefficients
       * @param point_index Index of reconstructed point
       * @param time_step_number Current time step number
       * @param cur_time Current time
       */
      void
      output_expansion_coefficients(
        const std::vector<double> &expansion_coefficients,
        const unsigned int         point_index,
        const unsigned int         time_step_number = 0,
        const double               cur_time         = 0.,
        const std::string         &function_symbol  = "f") const;
    };
  } // namespace VFP
} // namespace sapphirepp

#endif
