#ifndef VFP_PHASESPACERECONSTRUCTION_H
#define VFP_PHASESPACERECONSTRUCTION_H

#include <deal.II/lac/vector.h>

#include <array>
#include <filesystem>
#include <vector>

namespace sapphirepp
{
  namespace VFP
  {
    namespace PhaseSpaceReconstruction
    {
      std::vector<double>
      compute_phase_space_distribution(
        const std::vector<double>                      &mu_values,
        const std::vector<double>                      &phi_values,
        const std::vector<std::array<unsigned int, 3>> &lms_index_map,
        const dealii::Vector<double>                   &expansion_coefficients);

      void
      output_gnu_splot_data(const std::filesystem::path &path,
                            const std::vector<double>   &x_values,
                            const std::vector<double>   &y_values,
                            const std::vector<double>   &f_values);

      void
      output_gnu_splot_spherical_density_map(
        const std::filesystem::path &path,
        const std::vector<double>   &mu_values,
        const std::vector<double>   &phi_values,
        const std::vector<double>   &f_values);

      std::vector<double>
      create_range(const double       lower_bound,
                   const double       step_size,
                   const unsigned int n_intervals);
    } // namespace PhaseSpaceReconstruction
  }   // namespace VFP
      // namespace VFP
} // namespace sapphirepp

#endif
