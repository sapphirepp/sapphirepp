#include "phase-space-reconstruction.h"

#include <deal.II/lac/vector.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>



template <unsigned int dim>
sapphirepp::VFP::PhaseSpaceReconstruction<dim>::PhaseSpaceReconstruction(
  const VFPParameters<dim>                       &vfp_parameters,
  const std::vector<std::array<unsigned int, 3>> &lms_indices)
  : lms_indices{lms_indices}
  , theta_values{create_linear_range(0, M_PI, 76)}
  , phi_values{create_linear_range(0, 2 * M_PI, 76)}
  , real_spherical_harmonics{
      compute_real_spherical_harmonics(theta_values, phi_values, lms_indices)}
{
  static_cast<void>(vfp_parameters);
}



template <unsigned int dim>
std::vector<double>
sapphirepp::VFP::PhaseSpaceReconstruction<
  dim>::compute_phase_space_distribution(const dealii::Vector<double>
                                           &expansion_coefficients) const
{
  AssertDimension(expansion_coefficients.size(), lms_indices.size());

  std::vector<double> f(theta_values.size() * phi_values.size());

  for (std::size_t i = 0; i < lms_indices.size(); ++i)
    for (std::size_t j = 0; j < theta_values.size(); ++j)
      for (std::size_t k = 0; k < phi_values.size(); ++k)
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
std::vector<double>
sapphirepp::VFP::PhaseSpaceReconstruction<dim>::create_linear_range(
  const double       start,
  const double       stop,
  const unsigned int num)
{
  std::vector<double> values(num);
  for (unsigned int i = 0; i < num; ++i)
    values[i] = start + i * (stop - start) / (num - 1);
  // Sanitize
  values[0]       = start;
  values[num - 1] = stop;
  return values;
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
      const unsigned int m = lms_indices[i][1];
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