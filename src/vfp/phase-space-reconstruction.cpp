#include "phase-space-reconstruction.h"

#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>

#include "gsl/gsl_sf_legendre.h"

std::vector<double>
sapphirepp::VFP::PhaseSpace::compute_phase_space_distribution(
  const std::vector<double>                      &mu_values,
  const std::vector<double>                      &phi_values,
  const std::vector<std::array<unsigned int, 3>> &lms_index_map,
  const dealii::Vector<double>                   &expansion_coefficients)
{
  std::vector<double> cos_values(phi_values.size());
  std::vector<double> sin_values(phi_values.size());

  for (std::size_t i = 0; i < phi_values.size(); ++i)
    {
      cos_values[i] = std::cos(phi_values[i]);
      sin_values[i] = std::sin(phi_values[i]);
    }
  // NOTE: Type cast is necessary to choose a specific overload of std::cos
  // std::transform(phi_values.begin(), phi_values.end(), cos_values.begin(),
  //                static_cast<double (*)(double)>(std::cos));
  // std::transform(phi_values.begin(), phi_values.end(), sin_values.begin(),
  //                static_cast<double (*)(double)>(std::sin));

  std::vector<double> f(mu_values.size() * phi_values.size());
  for (std::size_t i = 0; i < expansion_coefficients.size(); ++i)
    {
      unsigned int l   = lms_index_map[i][0];
      unsigned int m   = lms_index_map[i][1];
      unsigned int s   = lms_index_map[i][2];
      auto         Plm = [l, m](auto &&PH1) {
        return gsl_sf_legendre_sphPlm(l, m, std::forward<decltype(PH1)>(PH1));
      };
      // auto         Plm = [l, m](double x) { gsl_sf_legendre_sphPlm(l, m, x);
      // };
      std::vector<double> Plm_values(mu_values.size());
      std::transform(mu_values.begin(),
                     mu_values.end(),
                     Plm_values.begin(),
                     Plm);

      for (std::size_t j = 0; j < mu_values.size(); ++j)
        {
          for (std::size_t k = 0; k < phi_values.size(); ++k)
            {
              if (m != 0)
                f[j * phi_values.size() + k] +=
                  expansion_coefficients[i] * //* (m % 2 ? -1. : 1.) * //
                                              // Condon-Shortley phase
                  Plm_values[j] * (s ? sin_values[k] : cos_values[k]);
              else
                f[j * phi_values.size() + k] +=
                  expansion_coefficients[i] // * (m % 2 ? -1. : 1.)
                  * Plm_values[j];
            }
        }
    }

  return f;
}


void
sapphirepp::VFP::PhaseSpace::output_gnu_splot_data(
  const std::filesystem::path &path,
  const std::vector<double>   &x_values,
  const std::vector<double>   &y_values,
  const std::vector<double>   &f_values)
{
  std::ofstream data_file(path);
  // See https://en.cppreference.com/w/cpp/types/numeric_limits/digits10
  data_file.precision(std::numeric_limits<double>::digits10);
  for (std::size_t i = 0; i < x_values.size(); ++i)
    {
      for (std::size_t j = 0; j < y_values.size(); ++j)
        data_file << x_values[i] << " " << y_values[j] << " "
                  << f_values[i * y_values.size() + j] << "\n";
      // Gnu plot data format requires the addition of an extra new line when
      // the x-coordinate changes. See
      // https://stackoverflow.com/questions/62729982/how-to-plot-a-3d-gnuplot-splot-surface-graph-with-data-from-a-file
      data_file << std::endl;
    }
}

void
sapphirepp::VFP::PhaseSpace::output_gnu_splot_spherical_density_map(
  const std::filesystem::path &path,
  const std::vector<double>   &mu_values,
  const std::vector<double>   &phi_values,
  const std::vector<double>   &f_values)
{
  std::ofstream data_file(path);
  data_file.precision(std::numeric_limits<double>::digits10);
  // x = mu
  double y = 0;
  double z = 0;
  for (std::size_t i = 0; i < mu_values.size(); ++i)
    {
      for (std::size_t j = 0; j < phi_values.size(); ++j)
        {
          y = std::sqrt(1 - mu_values[i] * mu_values[i]) *
              std::cos(phi_values[j]);
          z = std::sqrt(1 - mu_values[i] * mu_values[i]) *
              std::sin(phi_values[j]);
          data_file << mu_values[i] << " " << y << " " << z << " "
                    << f_values[i * phi_values.size() + j] << "\n";
        }
      data_file << std::endl;
    }
}

std::vector<double>
sapphirepp::VFP::PhaseSpace::create_range(const double       lower_bound,
                                          const double       step_size,
                                          const unsigned int n_intervals)
{
  std::vector<double> values(n_intervals + 1);
  double              start = lower_bound - step_size;
  std::generate(values.begin(), values.end(), [&start, step_size]() {
    start += step_size;
    return start;
  });
  return values;
}
