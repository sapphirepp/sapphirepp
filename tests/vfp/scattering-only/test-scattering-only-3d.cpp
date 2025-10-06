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
 * @file tests/vfp/scattering-only/test-scattering-only-3d.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement 3d tests for scattering-only example
 */

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <mpi.h>

#include "config.h"
#include "output-parameters.h"
#include "sapphirepp-logstream.h"
#include "test-run-vfp.h"
#include "vfp-parameters.h"
#include "vfp-solver.h"



const unsigned int dim = 3;



namespace sapinternal
{
  using namespace sapphirepp;

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
} // namespace sapinternal



int
main(int argc, char *argv[])
{
  try
    {
      /** [Main function setup] */
      using namespace sapphirepp;
      using namespace VFP;

      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);
      saplog.init(argc, argv);

      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      double max_L2_error = 0.;
      if (argc > 2)
        max_L2_error = std::stod(argv[2]);

      ParameterHandler        prm;
      PhysicalParameters      physical_parameters;
      Utils::OutputParameters output_parameters;
      VFPParameters<dim>      vfp_parameters(vfp_flags);

      physical_parameters.declare_parameters(prm);
      output_parameters.declare_parameters(prm);
      vfp_parameters.declare_parameters(prm);

      prm.parse_input(parameter_filename);

      physical_parameters.parse_parameters(prm);
      output_parameters.parse_parameters(prm);
      vfp_parameters.parse_parameters(prm);
      /** [Main function setup] */


      /** [Run simulation] */
      VFPSolver<dim> vfp_solver(vfp_parameters,
                                physical_parameters,
                                output_parameters);
      vfp_solver.run();
      /** [Run simulation] */


      /** [Compare to exact solution] */
      saplog << "Compare to exact solution" << std::endl;

      sapinternal::AnalyticSolution<dim> exact_solution(
        physical_parameters,
        vfp_solver.get_pde_system().system_size,
        vfp_parameters.final_time);

      test_run_vfp_output<dim>(vfp_solver,
                               vfp_parameters,
                               output_parameters,
                               exact_solution,
                               0,
                               vfp_parameters.final_time,
                               "exact_solution");


      test_run_vfp_error<dim>(vfp_solver, exact_solution, saplog, max_L2_error);

      sapphirepp::saplog << "Succeeded test run VFP." << std::endl;
      /** [Compare to exact solution] */
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
