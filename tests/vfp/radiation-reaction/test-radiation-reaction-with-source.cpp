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
 * @file
 * tests/vfp/radiation-reaction/test-radiation-reaction-with-rotation.cpp
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @brief Implement tests for radiation-reaction example with source
 */

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <mpi.h>

#include <cmath>

#include "config.h"
#include "output-parameters.h"
#include "sapphirepp-logstream.h"
#include "test-run-vfp.h"
#include "vfp-parameters.h"
#include "vfp-solver.h"



const unsigned int dim = sapphirepp::VFP::dimension;



namespace sapinternal
{
  using namespace dealii;
  using namespace sapphirepp;

  template <unsigned int dim>
  class AnalyticSolution : public dealii::Function<dim>
  {
  public:
    AnalyticSolution(const PhysicalParameters      &physical_parameters,
                     const VFP::VFPParameters<dim> &vfp_parameters,
                     const unsigned int             system_size)
      : Function<dim>(system_size)
      , prm{physical_parameters}
      , vfp_parameters{vfp_parameters}
      , lms_indices{VFP::PDESystem::create_lms_indices(system_size)}
    {}



    void
    vector_value(const Point<dim> &point, Vector<double> &g) const override
    {
      AssertDimension(g.size(), this->n_components);

      g = 0;
      // Momentum coordinate
      const double p = std::exp(point[0]);

      // Time
      const double t = this->get_time();
      const double tau_R =
        vfp_parameters.reference_units.radiation_reaction_characteristic_time /
        (std::pow(vfp_parameters.mass, -3) *
         std::pow(vfp_parameters.charge, 4));


      const double tau_inj = tau_R / (prm.B0 * prm.B0 * prm.p_inj);
      const double tau_c   = tau_R / (prm.B0 * prm.B0 * p);

      // Analytical solution
      // std::cout << (t >= (tau_c - tau_inj))  << "\n";
      // t >= (tau_c - tau_inj)
      if (t >= (tau_c - tau_inj) and p <= prm.p_inj)
        g[0] = prm.injection_rate / (std::sqrt(4 * std::numbers::pi) * p) *
               prm.p_inj * tau_inj;
    }


  private:
    const PhysicalParameters                       prm;
    const VFP::VFPParameters<dim>                  vfp_parameters;
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

      dealii::ParameterHandler prm;
      PhysicalParameters       physical_parameters;
      Utils::OutputParameters  output_parameters;
      VFPParameters<dim>       vfp_parameters(vfp_flags);

      physical_parameters.declare_parameters(prm);
      output_parameters.declare_parameters(prm);
      vfp_parameters.declare_parameters(prm);

      prm.parse_input(parameter_filename);

      physical_parameters.parse_parameters(prm);
      output_parameters.parse_parameters(prm);
      vfp_parameters.parse_parameters(prm);
      /** [Main function setup] */


      /** [Setup exact solution] */
      const unsigned int system_size = (vfp_parameters.expansion_order + 1) *
                                       (vfp_parameters.expansion_order + 1);
      sapinternal::AnalyticSolution<dim> exact_solution(physical_parameters,
                                                        vfp_parameters,
                                                        system_size);
      /** [Setup exact solution] */

      /** [Start test run] */
      return test_run_vfp<dim>(vfp_parameters,
                               physical_parameters,
                               output_parameters,
                               exact_solution,
                               max_L2_error,
                               nullptr,
                               false);
      /** [Start test run] */
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
                << "Unknown exception!"
                << "\n"
                << "Aborting!"
                << "\n"
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
