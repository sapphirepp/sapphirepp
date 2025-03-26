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
 * @file probe-location.cpp
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de),
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement @ref sapphirepp::VFP::ProbeLocation
 */

#include "probe-location.h"

#include <deal.II/base/array_view.h>
#include <deal.II/base/mpi_remote_point_evaluation.h>

#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools_evaluate.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <cmath>
#include <fstream>
#include <limits>

#include "tools.h"



template <unsigned int dim>
sapphirepp::VFP::ProbeLocation<dim>::ProbeLocation(
  const VFPParameters<dim>                       &vfp_parameters,
  const Utils::OutputParameters                  &output_parameters,
  const std::vector<std::array<unsigned int, 3>> &lms_indices)
  : output_parameters{output_parameters}
  , lms_indices{lms_indices}
  , system_size{static_cast<unsigned int>(lms_indices.size())}
  , perform_probe_location{vfp_parameters.perform_probe_location}
  , perform_phase_space_reconstruction{vfp_parameters
                                         .perform_phase_space_reconstruction}
  , theta_values{Utils::Tools::create_linear_range(0,
                                                   M_PI,
                                                   vfp_parameters.n_theta)}
  , phi_values{Utils::Tools::create_linear_range(0,
                                                 2 * M_PI,
                                                 vfp_parameters.n_phi)}
  , real_spherical_harmonics{compute_real_spherical_harmonics(theta_values,
                                                              phi_values,
                                                              lms_indices)}
  , rpe_cache(1e-6, true)
{
  LogStream::Prefix p0("VFP", saplog);
  LogStream::Prefix p("ProbeLocation", saplog);
  LogStream::Prefix p1("Constructor", saplog);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    probe_location_points = vfp_parameters.probe_location_points;

  if (probe_location_points.size() > 0)
    {
      saplog << "Probe location will be performed at the following points:"
             << std::endl;
      for (const auto &point : probe_location_points)
        saplog << "  " << point << std::endl;
    }
}



template <unsigned int dim>
void
sapphirepp::VFP::ProbeLocation<dim>::reinit(
  const Triangulation<dim> &triangulation,
  const Mapping<dim>       &mapping)
{
  if (!perform_probe_location)
    return;
  LogStream::Prefix p("ProbeLocation", saplog);
  saplog << "Setting up points for probe location" << std::endl;

  rpe_cache.reinit(probe_location_points, triangulation, mapping);
  AssertThrow(rpe_cache.all_points_found(),
              ExcMessage("Point for reconstruction is outside domain."));
  Assert(rpe_cache.is_map_unique(),
         ExcMessage("Mapping between points and cells in not unique."));
}



template <unsigned int dim>
void
sapphirepp::VFP::ProbeLocation<dim>::probe_all_points(
  const DoFHandler<dim>            &dof_handler,
  const Mapping<dim>               &mapping,
  const PETScWrappers::MPI::Vector &solution,
  const unsigned int                time_step_number,
  const double                      cur_time) const
{
  if (!perform_probe_location)
    return;
  LogStream::Prefix p("ProbeLocation", saplog);
  saplog << "Probe points in reduced phase space" << std::endl;
  LogStream::Prefix p2("Evaluate", saplog);
  Assert(rpe_cache.is_ready(), ExcMessage("rpe_cache is no initialized."));

  using CellData =
    typename dealii::Utilities::MPI::RemotePointEvaluation<dim, dim>::CellData;


  const std::function<void(const ArrayView<std::vector<double>> &,
                           const CellData &)>
    evaluate_function = [&](const ArrayView<std::vector<double>> &coefficients,
                            const CellData                       &cell_data) {
      std::vector<Point<dim>>              unit_points;
      std::vector<types::global_dof_index> local_dof_indices_vector;

      for (const auto cell_index : cell_data.cell_indices())
        {
          const auto &cell = cell_data.get_active_cell_iterator(cell_index)
                               ->as_dof_handler_iterator(dof_handler);

          // Convert unit point array view to vector
          unit_points.clear();
          for (const Point<dim> &p : cell_data.get_unit_points(cell_index))
            unit_points.push_back(p);

          // Initialize FEValues for point evaluation
          FEValues<dim> fe_values(mapping,
                                  cell->get_fe(),
                                  Quadrature<dim>(unit_points),
                                  update_values);
          fe_values.reinit(cell);

          // Resize coefficient buffer
          auto local_coefficients =
            cell_data.get_data_view(cell_index, coefficients);
          for (auto &coeff : local_coefficients)
            coeff.resize(system_size);

          // Get local_dof_indices as ArrayView
          local_dof_indices_vector.resize(cell->get_fe().n_dofs_per_cell());
          cell->get_dof_indices(local_dof_indices_vector);
          const auto local_dof_indices =
            make_array_view(local_dof_indices_vector);

          fe_values.get_function_values(solution,
                                        local_dof_indices,
                                        local_coefficients,
                                        false);
        }
    };


  const std::vector<std::vector<double>> coefficients_at_all_points =
    rpe_cache.evaluate_and_process(evaluate_function);
  AssertDimension(coefficients_at_all_points.size(),
                  probe_location_points.size());

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      // loop over all points and compute the phase reconstruction
      unsigned int point_counter = 0;
      for (const auto &expansion_coefficients : coefficients_at_all_points)
        {
          output_f_lms(expansion_coefficients,
                       point_counter,
                       time_step_number,
                       cur_time);

          if (perform_phase_space_reconstruction)
            {
              std::vector<double> f_values =
                compute_phase_space_distribution(expansion_coefficients);

              output_gnu_splot_data(f_values,
                                    point_counter,
                                    time_step_number,
                                    cur_time);
              output_gnu_splot_spherical_density_map(f_values,
                                                     point_counter,
                                                     time_step_number,
                                                     cur_time);
            }

          ++point_counter;
        }
    }
}



template <unsigned int dim>
std::vector<double>
sapphirepp::VFP::ProbeLocation<dim>::compute_phase_space_distribution(
  const std::vector<double> &expansion_coefficients) const
{
  AssertDimension(expansion_coefficients.size(), lms_indices.size());

  std::vector<double> f(theta_values.size() * phi_values.size());

  for (unsigned int i = 0; i < lms_indices.size(); ++i)
    for (unsigned int j = 0; j < theta_values.size(); ++j)
      for (unsigned int k = 0; k < phi_values.size(); ++k)
        f[j * phi_values.size() + k] +=
          expansion_coefficients[i] * real_spherical_harmonics(j, k, i);

  return f;
}



template <unsigned int dim>
void
sapphirepp::VFP::ProbeLocation<dim>::output_gnu_splot_data(
  const std::vector<double> &f_values,
  const unsigned int         point_index,
  const unsigned int         time_step_number,
  const double               cur_time) const
{
  AssertIndexRange(point_index, probe_location_points.size());

  const std::string filename =
    "surface_plot_distribution_function_point_" +
    Utilities::int_to_string(point_index, 2) + "_t_" +
    Utilities::int_to_string(time_step_number,
                             output_parameters.n_digits_for_counter) +
    ".dat";

  saplog << "Write reconstruction at (" << probe_location_points[point_index]
         << ") to file " << filename << std::endl;
  std::ofstream data_file(output_parameters.output_path / filename);
  // See https://en.cppreference.com/w/cpp/types/numeric_limits/digits10
  data_file.precision(std::numeric_limits<double>::digits10);

  data_file << "# f(t, x, p, theta, mu) at (x,|p|) = ("
            << probe_location_points[point_index] << ") for t = " << cur_time
            << "\n";
  data_file << "# theta phi f(theta, phi) \n";
  for (unsigned int i = 0; i < theta_values.size(); ++i)
    {
      for (unsigned int j = 0; j < phi_values.size(); ++j)
        data_file << theta_values[i] << " " << phi_values[j] << " "
                  << f_values[i * theta_values.size() + j] << "\n";
      // Gnu plot data format requires the addition of an extra new line when
      // the x-coordinate changes. See
      // https://stackoverflow.com/questions/62729982/how-to-plot-a-3d-gnuplot-splot-surface-graph-with-data-from-a-file
      data_file << std::endl;
    }
}



template <unsigned int dim>
void
sapphirepp::VFP::ProbeLocation<dim>::output_gnu_splot_spherical_density_map(
  const std::vector<double> &f_values,
  const unsigned int         point_index,
  const unsigned int         time_step_number,
  const double               cur_time) const
{
  AssertIndexRange(point_index, probe_location_points.size());

  const std::string filename =
    "spherical_density_map_point_" + Utilities::int_to_string(point_index, 2) +
    "_t_" +
    Utilities::int_to_string(time_step_number,
                             output_parameters.n_digits_for_counter) +
    ".dat";

  saplog << "Write reconstruction at (" << probe_location_points[point_index]
         << ") to file " << filename << std::endl;
  std::ofstream data_file(output_parameters.output_path / filename);
  data_file.precision(std::numeric_limits<double>::digits10);

  data_file << "# f(t, x, p) at (x,|p|) = ("
            << probe_location_points[point_index] << ") for t = " << cur_time
            << "\n";
  data_file << "# n_{p_x} n_{p_y} n_{p_z} f(p_x, p_z, p_y) \n";
  for (unsigned int i = 0; i < theta_values.size(); ++i)
    {
      for (unsigned int j = 0; j < phi_values.size(); ++j)
        {
          const double x = std::cos(theta_values[i]);
          const double y = std::sin(theta_values[i]) * std::cos(phi_values[j]);
          const double z = std::sin(theta_values[i]) * std::sin(phi_values[j]);
          data_file << x << " " << y << " " << z << " "
                    << f_values[i * phi_values.size() + j] << "\n";
        }
      data_file << std::endl;
    }
}

template <unsigned int dim>
void
sapphirepp::VFP::ProbeLocation<dim>::output_f_lms(
  const std::vector<double> &expansion_coefficients,
  const unsigned int         point_index,
  const unsigned int         time_step_number,
  const double               cur_time) const
{
  AssertIndexRange(point_index, probe_location_points.size());

  const std::string filename = "f_lms_values_at_point_" +
                               Utilities::int_to_string(point_index, 2) +
                               ".dat";

  saplog << "Write values of expansion coefficients at point ("
         << probe_location_points[point_index] << ") to file " << filename
         << std::endl;
  std::ofstream data_file(output_parameters.output_path / filename,
                          time_step_number == 0 ? std::ios_base::out :
                                                  std::ios_base::app);

  // See https://en.cppreference.com/w/cpp/types/numeric_limits/digits10
  std::stringstream sstream;
  sstream.precision(std::numeric_limits<double>::digits10);

  if (time_step_number == 0)
    {
      sstream << "# f(t, x, p ) at (x,|p|) = ("
              << probe_location_points[point_index] << ")"
              << "\n";
      sstream << "# time_step_number cur_time ";
      for (auto &lms : lms_indices)
        sstream << "f_" << lms[0] << lms[1] << lms[2] << " ";
      sstream << "\n";
    }
  sstream << time_step_number << " " << cur_time << " ";
  for (auto f_lms : expansion_coefficients)
    sstream << f_lms << " ";
  // Write string to  data file

  data_file << sstream.str() << "\n";
}



template <unsigned int dim>
dealii::Table<3, double>
sapphirepp::VFP::ProbeLocation<dim>::compute_real_spherical_harmonics(
  const std::vector<double>                      &theta_values,
  const std::vector<double>                      &phi_values,
  const std::vector<std::array<unsigned int, 3>> &lms_indices)
{
  dealii::Table<3, double> y_lms(theta_values.size(),
                                 phi_values.size(),
                                 lms_indices.size());

  for (unsigned int i = 0; i < lms_indices.size(); ++i)
    {
      const unsigned int l = lms_indices[i][0];
      const int          m = static_cast<int>(lms_indices[i][1]);
      const unsigned int s = lms_indices[i][2];

      for (unsigned int j = 0; j < theta_values.size(); ++j)
        {
          for (unsigned int k = 0; k < phi_values.size(); ++k)
            {
              switch (s)
                {
                  case 0:
                    {
                      y_lms(j, k, i) = boost::math::spherical_harmonic_r(
                        l, m, theta_values[j], phi_values[k]);
                      if (m > 0)
                        y_lms(j, k, i) *= M_SQRT2;
                      break;
                    }
                  case 1:
                    {
                      Assert(m > 0, ExcInvalidState());
                      y_lms(j, k, i) =
                        M_SQRT2 * boost::math::spherical_harmonic_i(
                                    l, m, theta_values[j], phi_values[k]);
                      break;
                    }
                  default:
                    Assert(false, ExcInvalidState());
                }
            }
        }
    }

  return y_lms;
}



// Explicit instantiation
template class sapphirepp::VFP::ProbeLocation<1>;
template class sapphirepp::VFP::ProbeLocation<2>;
template class sapphirepp::VFP::ProbeLocation<3>;
