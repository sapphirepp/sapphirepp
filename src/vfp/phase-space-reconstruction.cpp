#include "phase-space-reconstruction.h"

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools_evaluate.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>

#include "tools.h"



template <unsigned int dim>
sapphirepp::VFP::PhaseSpaceReconstruction<dim>::PhaseSpaceReconstruction(
  const VFPParameters<dim>                       &vfp_parameters,
  const Utils::OutputParameters                  &output_parameters,
  const std::vector<std::array<unsigned int, 3>> &lms_indices)
  : output_parameters{output_parameters}
  , lms_indices{lms_indices}
  , perform_phase_space_reconstruction{vfp_parameters
                                         .perform_phase_space_reconstruction}
  , theta_values{Utils::Tools::create_linear_range(0,
                                                   M_PI,
                                                   vfp_parameters.n_theta)}
  , phi_values{Utils::Tools::create_linear_range(0,
                                                 2 * M_PI,
                                                 vfp_parameters.n_phi)}
  , real_spherical_harmonics{
      compute_real_spherical_harmonics(theta_values, phi_values, lms_indices)}
{
  LogStream::Prefix p("PSRec", saplog);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    reconstruction_points = vfp_parameters.reconstruction_points;

  if (reconstruction_points.size() > 0)
    {
      saplog << "Phase space will be reconstructed at:" << std::endl;
      for (const auto &point : reconstruction_points)
        saplog << "  " << point << std::endl;
      saplog << std::endl;
    }
}



template <unsigned int dim>
void
sapphirepp::VFP::PhaseSpaceReconstruction<dim>::reinit(
  const Triangulation<dim> &triangulation,
  const Mapping<dim>       &mapping)
{
  if (!perform_phase_space_reconstruction)
    return;
  LogStream::Prefix p("PSRec", saplog);
  saplog << "Reconstruction of phase space is set up" << std::endl;

  rpe_cache.reinit(reconstruction_points, triangulation, mapping);
}



template <unsigned int dim>
void
sapphirepp::VFP::PhaseSpaceReconstruction<dim>::reconstruct_all_points(
  const DoFHandler<dim>            &dof_handler,
  const PETScWrappers::MPI::Vector &solution,
  const unsigned int                time_step_number) const
{
  if (!perform_phase_space_reconstruction)
    return;
  LogStream::Prefix p("PSRec", saplog);
  saplog << "Perform reconstruction of phase space" << std::endl;

  // NOTE: point_values needs the template argument n_components, which needs to
  // be known at compile time. Unfortunately it is not possible to compute at
  // compile time, because l_max is a run-time parameter.
  //
  // TODO: Find a workaround. Probably a lambda that is passed to
  // rpe.evaluate_and_process() that does not use FEPointEvaluation, which
  // requires the template argument.
  constexpr unsigned int n_components = (3 + 1) * (3 + 1);
  // Return value: std::vector<dealii::Tensor<1,n_components,double>>
  const auto coefficients_at_all_points =
    VectorTools::point_values<n_components>(rpe_cache, dof_handler, solution);

  // NOTE: In an actual implementaion, the computation of phi and mu ranges and
  // also the values of cos phi, sin phi and Plm(mu) should only be done once.
  // For example, in the setup up phase of a PostProcessing Object.
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      // loop over all points and compute the phase reconstruction
      unsigned int point_counter = 0;
      for (const auto &expansion_coefficients : coefficients_at_all_points)
        {
          // TODO: In the actual implemenation compute phase distribution should
          // directly work with the tensor. There is no need to unroll.
          dealii::Vector<double> flms_values(n_components);
          expansion_coefficients.unroll(flms_values.begin(), flms_values.end());
          std::vector<double> f_values =
            compute_phase_space_distribution(flms_values);

          std::string path = output_parameters.output_path;
          output_gnu_splot_data(
            path + "/surface_plot_distribution_function_p_" +
              Utilities::to_string(point_counter) + "_at_t_" +
              Utilities::to_string(time_step_number) + ".dat",
            f_values);
          output_gnu_splot_spherical_density_map(
            path + "/spherical_density_map_p_" +
              Utilities::to_string(point_counter) + "_at_t_" +
              Utilities::to_string(time_step_number) + ".dat",
            f_values);
          ++point_counter;
        }
    }
}



template <unsigned int dim>
std::vector<double>
sapphirepp::VFP::PhaseSpaceReconstruction<
  dim>::compute_phase_space_distribution(const dealii::Vector<double>
                                           &expansion_coefficients) const
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
sapphirepp::VFP::PhaseSpaceReconstruction<dim>::output_gnu_splot_data(
  const std::filesystem::path &path,
  const std::vector<double>   &f_values) const
{
  std::ofstream data_file(path);
  // See https://en.cppreference.com/w/cpp/types/numeric_limits/digits10
  data_file.precision(std::numeric_limits<double>::digits10);
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
sapphirepp::VFP::PhaseSpaceReconstruction<dim>::
  output_gnu_splot_spherical_density_map(
    const std::filesystem::path &path,
    const std::vector<double>   &f_values) const
{
  std::ofstream data_file(path);
  data_file.precision(std::numeric_limits<double>::digits10);
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
dealii::Table<3, double>
sapphirepp::VFP::PhaseSpaceReconstruction<dim>::
  compute_real_spherical_harmonics(
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
                    Assert(false, ExcInvalidState())
                }
            }
        }
    }

  return y_lms;
}



// Explicit instantiation
template class sapphirepp::VFP::PhaseSpaceReconstruction<1>;
template class sapphirepp::VFP::PhaseSpaceReconstruction<2>;
template class sapphirepp::VFP::PhaseSpaceReconstruction<3>;