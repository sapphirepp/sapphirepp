#ifndef VFP_PHASESPACERECONSTRUCTION_H
#define VFP_PHASESPACERECONSTRUCTION_H

#include <deal.II/base/table.h>

#include <deal.II/lac/vector.h>

#include <array>
#include <filesystem>
#include <vector>

#include "vfp-parameters.h"

namespace sapphirepp
{
  namespace VFP
  {
    template <unsigned int dim>
    class PhaseSpaceReconstruction
    {
    public:
      PhaseSpaceReconstruction(
        const VFPParameters<dim>                       &vfp_parameters,
        const std::vector<std::array<unsigned int, 3>> &lms_indices);



      std::vector<double>
      compute_phase_space_distribution(
        const dealii::Vector<double> &expansion_coefficients) const;



      void
      output_gnu_splot_data(const std::filesystem::path &path,
                            const std::vector<double>   &f_values) const;



      void
      output_gnu_splot_spherical_density_map(
        const std::filesystem::path &path,
        const std::vector<double>   &f_values) const;



      static std::vector<double>
      create_linear_range(const double       start,
                          const double       stop,
                          const unsigned int num);



      static dealii::Table<3, double>
      compute_real_spherical_harmonics(
        const std::vector<double>                      &theta_values,
        const std::vector<double>                      &phi_values,
        const std::vector<std::array<unsigned int, 3>> &lms_indices);



    private:
      const std::vector<std::array<unsigned int, 3>> lms_indices;

      const std::vector<double> theta_values;
      const std::vector<double> phi_values;

      const dealii::Table<3, double> real_spherical_harmonics;
    };
  } // namespace VFP
} // namespace sapphirepp

#endif
