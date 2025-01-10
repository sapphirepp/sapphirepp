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
 * @file examples/scattering-only/scattering-only.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement main function for scattering-only example
 */

/** [Includes] */
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <mpi.h>

#include "config.h"
#include "output-parameters.h"
#include "sapphirepp-logstream.h"
#include "vfp-parameters.h"
#include "vfp-solver.h"
/** [Includes] */


/** [AnalyticSolution namespace] */
namespace sapinternal
{
  namespace AnalyticSolutionImplementation
  {
    using namespace sapphirepp;
    /** [AnalyticSolution namespace] */

    /** [AnalyticSolution constructor] */
    template <unsigned int dim>
    class AnalyticSolution : public dealii::Function<dim>
    {
    public:
      AnalyticSolution(const PhysicalParameters &physical_parameters,
                       const unsigned int        system_size,
                       const double              time)
        : dealii::Function<dim>(system_size, time)
        , prm{physical_parameters}
        , lms_indices{VFP::PDESystem::create_lms_indices(system_size)}
      {}
      /** [AnalyticSolution constructor] */



      /** [AnalyticSolution value] */
      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        for (unsigned int i = 0; i < f.size(); ++i)
          {
            const unsigned int l = lms_indices[i][0];
            const double       t = this->get_time();

            f[i] = prm.f0 * std::exp(-prm.nu * l * (l + 1) / 2. * t);
          }
      }



    private:
      const PhysicalParameters                       prm;
      const std::vector<std::array<unsigned int, 3>> lms_indices;
    };
  } // namespace AnalyticSolutionImplementation
} // namespace sapinternal
/** [AnalyticSolution value] */



/** [Main function] */
int
main(int argc, char *argv[])
{
  /** [Main function] */
  /** [Try-Catch begin] */
  try
    {
      using namespace sapphirepp;
      using namespace VFP;
      /** [Try-Catch begin] */
      /** [MPI initialization] */
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);
      /** [MPI initialization] */
      /** [Saplog] */

      saplog.init(argc, argv);
      /** [Saplog] */

      /** [Command line argument] */
      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];
      /** [Command line argument] */

      /** [Run time parameters] */
      dealii::ParameterHandler prm;
      PhysicalParameters       physical_parameters;
      Utils::OutputParameters  output_parameters;
      VFPParameters<dimension> vfp_parameters;
      /** [Run time parameters] */

      /** [Declare parameters] */
      physical_parameters.declare_parameters(prm);
      output_parameters.declare_parameters(prm);
      vfp_parameters.declare_parameters(prm);
      /** [Declare parameters] */

      /** [Parse parameters] */
      prm.parse_input(parameter_filename);

      physical_parameters.parse_parameters(prm);
      output_parameters.parse_parameters(prm);
      vfp_parameters.parse_parameters(prm);
      /** [Parse parameters] */

      /** [VFP Solver] */
      VFPSolver<dimension> vfp_solver(vfp_parameters,
                                      physical_parameters,
                                      output_parameters);
      vfp_solver.run();
      /** [VFP Solver] */

      /** [Create AnalyticSolution] */
      using namespace sapinternal::AnalyticSolutionImplementation;

      AnalyticSolution<dimension> analytic_solution(
        physical_parameters,
        vfp_solver.get_pde_system().system_size,
        vfp_parameters.final_time);
      /** [Create AnalyticSolution] */

      /** [Calculate L2-error] */
      const double L2_error =
        vfp_solver.compute_global_error(analytic_solution,
                                        dealii::VectorTools::L2_norm,
                                        dealii::VectorTools::L2_norm);
      const double L2_norm =
        vfp_solver.compute_weighted_norm(dealii::VectorTools::L2_norm,
                                         dealii::VectorTools::L2_norm);

      saplog << "L2_error = " << L2_error << ", L2_norm = " << L2_norm
             << ", rel error = " << L2_error / L2_norm << std::endl;
      /** [Calculate L2-error] */

      /** [Calculate analytic solution] */
      PETScWrappers::MPI::Vector analytic_solution_vector;
      analytic_solution_vector.reinit(
        vfp_solver.get_dof_handler().locally_owned_dofs(), MPI_COMM_WORLD);
      dealii::VectorTools::interpolate(vfp_solver.get_dof_handler(),
                                       analytic_solution,
                                       analytic_solution_vector);
      /** [Calculate analytic solution] */

      /** [Output analytic solution] */
      dealii::DataOut<dimension> data_out;
      data_out.attach_dof_handler(vfp_solver.get_dof_handler());
      data_out.add_data_vector(analytic_solution_vector,
                               PDESystem::create_component_name_list(
                                 vfp_solver.get_pde_system().system_size,
                                 "analytic_f_"));
      data_out.build_patches(vfp_parameters.polynomial_degree);
      output_parameters.write_results<dimension, dimension>(
        data_out, 0, vfp_parameters.final_time, "analytic_solution");
      /** [Output analytic solution] */
      /** [Try-Catch end] */
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
/** [Try-Catch end] */
