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
 * @file tests/vfpscattering-only/test-scattering-only-1d.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement 1d tests for scattering-only example
 */

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <mpi.h>

#include <fstream>

#include "config.h"
#include "output-parameters.h"
#include "sapphirepp-logstream.h"
#include "vfp-parameters.h"
#include "vfp-solver.h"


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



      void
      vector_value([[maybe_unused]] const dealii::Point<dim> &point,
                   dealii::Vector<double>                    &f) const override
      {
        AssertDimension(f.size(), this->n_components);

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


const unsigned int dim = 1;

int
convergence_with_expansion_order(const std::string         &parameter_filename,
                                 const std::vector<double> &values)
{
  using namespace sapphirepp;
  using namespace VFP;

  saplog.push("Tests");
  saplog << "Compute convergence with expansion_order" << std::endl;

  Timer                   timer;
  ParameterHandler        prm;
  VFPParameters<dim>      vfp_parameters(vfp_flags);
  PhysicalParameters      physical_parameters;
  Utils::OutputParameters output_parameters;

  vfp_parameters.declare_parameters(prm);
  physical_parameters.declare_parameters(prm);
  output_parameters.declare_parameters(prm);

  prm.parse_input(parameter_filename);

  vfp_parameters.parse_parameters(prm);
  physical_parameters.parse_parameters(prm);
  output_parameters.parse_parameters(prm);

  std::ofstream log_file(output_parameters.results_path /
                           output_parameters.simulation_id /
                           "convergence_expansion_order.csv",
                         std::ios::app);
  saplog.attach(log_file, false);

  saplog.pop();
  saplog << "# "
         << "expansion_order"
         << "\t L2\t Linfty\t CPU time [s]\t Wall time [s]\t n_dof"
         << std::endl;
  saplog.push("Tests");
  for (const double value : values)
    {
      saplog << "expansion_order"
             << "=" << value << std::endl;
      vfp_parameters.expansion_order = static_cast<unsigned int>(value);
      output_parameters.base_file_name =
        "expansion_order_" + dealii::Utilities::to_string(value);

      VFPSolver<dim> vfp_solver(vfp_parameters,
                                physical_parameters,
                                output_parameters);
      vfp_solver.run();
      timer.stop();

      using namespace sapinternal::AnalyticSolutionImplementation;
      AnalyticSolution<dimension> analytic_solution(
        physical_parameters,
        vfp_solver.get_pde_system().system_size,
        vfp_parameters.final_time);

      const double L2_error =
        vfp_solver.compute_global_error(analytic_solution,
                                        dealii::VectorTools::L2_norm,
                                        dealii::VectorTools::L2_norm);
      const double L2_norm =
        vfp_solver.compute_weighted_norm(dealii::VectorTools::L2_norm,
                                         dealii::VectorTools::L2_norm);

      saplog << "L2_error = " << L2_error << ", L2_norm = " << L2_norm
             << ", rel error = " << L2_error / L2_norm
             << ", CPU/wall time = " << timer.cpu_time() << "/"
             << timer.wall_time() << " s" << std::endl;
      saplog.push("Tests");
    }

  saplog.detach();
  log_file.close();

  return 0;
}


int
test_run(const std::string &parameter_filename, const double max_L2_error)
{
  using namespace sapphirepp;
  using namespace VFP;

  saplog.push("Tests");
  saplog << "Test run with parameter file \"" << parameter_filename
         << "\" and maximal L2 error of " << max_L2_error << std::endl;

  dealii::Timer           timer;
  ParameterHandler        prm;
  VFPParameters<dim>      vfp_parameters(vfp_flags);
  PhysicalParameters      physical_parameters;
  Utils::OutputParameters output_parameters;

  vfp_parameters.declare_parameters(prm);
  physical_parameters.declare_parameters(prm);
  output_parameters.declare_parameters(prm);
  prm.parse_input(parameter_filename);

  vfp_parameters.parse_parameters(prm);
  physical_parameters.parse_parameters(prm);
  output_parameters.parse_parameters(prm);

  timer.start();
  VFPSolver<dim> vfp_solver(vfp_parameters,
                            physical_parameters,
                            output_parameters);
  vfp_solver.run();
  timer.stop();

  using namespace sapinternal::AnalyticSolutionImplementation;
  AnalyticSolution<dimension> analytic_solution(
    physical_parameters,
    vfp_solver.get_pde_system().system_size,
    vfp_parameters.final_time);

  const double L2_error =
    vfp_solver.compute_global_error(analytic_solution,
                                    dealii::VectorTools::L2_norm,
                                    dealii::VectorTools::L2_norm);
  const double L2_norm =
    vfp_solver.compute_weighted_norm(dealii::VectorTools::L2_norm,
                                     dealii::VectorTools::L2_norm);

  saplog << "L2_error = " << L2_error << ", L2_norm = " << L2_norm
         << ", rel error = " << L2_error / L2_norm
         << ", CPU/wall time = " << timer.cpu_time() << "/" << timer.wall_time()
         << " s" << std::endl;

  AssertThrow(L2_error / L2_norm < max_L2_error,
              dealii::ExcMessage(
                "L2 error is too large! (" +
                dealii::Utilities::to_string(L2_error / L2_norm) + " > " +
                dealii::Utilities::to_string(max_L2_error) + ")"));
  saplog.pop();
  return 0;
}


int
main(int argc, char *argv[])
{
  try
    {
      using namespace sapphirepp;
      using namespace VFP;
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);

      saplog.init();
      saplog.pop();
      saplog.depth_console(2);
      saplog.depth_file(0);
      saplog << "Start test-scattering-only-1d" << std::endl;

      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      std::string run_type = "test_run";
      if (argc > 2)
        run_type = argv[2];

      if (run_type == "test_run")
        {
          double max_L2_error = VFPParameters<dimension>::epsilon_d;
          if (argc > 3)
            max_L2_error = std::stod(argv[3]);
          test_run(parameter_filename, max_L2_error);
        }
      else if (run_type == "expansion_order")
        {
          std::vector<double> values = {1, 2, 3, 4, 5, 6, 7, 8};
          if (argc > 3)
            {
              values.clear();
              for (int i = 3; i < argc; i++)
                values.push_back(std::stod(argv[i]));
            }
          convergence_with_expansion_order(parameter_filename, values);
        }
      else
        AssertThrow(false,
                    dealii::ExcMessage("Unknown run type \"" + run_type +
                                       "\"!"));
    }
  catch (std::exception &exc)
    {
      sapphirepp::saplog.print_error(exc);
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl;
      std::cerr << "\n"
                << "----------------------------------------------------"
                << "\n"
                << "Unknown exception!" << "\n"
                << "Aborting!" << "\n"
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
