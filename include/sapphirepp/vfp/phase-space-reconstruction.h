#ifndef VFP_PHASESPACERECONSTRUCTION_H
#define VFP_PHASESPACERECONSTRUCTION_H

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

    template <unsigned int dim>
    class PhaseSpaceReconstruction
    {
    public:
      PhaseSpaceReconstruction(
        const VFPParameters<dim>                       &vfp_parameters,
        const Utils::OutputParameters                  &output_parameters,
        const std::vector<std::array<unsigned int, 3>> &lms_indices);



      void
      reinit(const Triangulation<dim> &triangulation,
             const Mapping<dim>       &mapping);



      void
      reconstruct_all_points(const DoFHandler<dim>            &dof_handler,
                             const PETScWrappers::MPI::Vector &solution,
                             const unsigned int time_step_number = 0) const;



      std::vector<double>
      compute_phase_space_distribution(
        const Vector<double> &expansion_coefficients) const;



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



      static Table<3, double>
      compute_real_spherical_harmonics(
        const std::vector<double>                      &theta_values,
        const std::vector<double>                      &phi_values,
        const std::vector<std::array<unsigned int, 3>> &lms_indices);



    private:
      /** Output parameter */
      const Utils::OutputParameters output_parameters;

      const std::vector<std::array<unsigned int, 3>> lms_indices;

      const std::vector<double> theta_values;
      const std::vector<double> phi_values;

      const Table<3, double> real_spherical_harmonics;

      std::vector<Point<dim>> reconstruction_points;

      Utilities::MPI::RemotePointEvaluation<dim, dim> rpe_cache;
    };
  } // namespace VFP
} // namespace sapphirepp

#endif
