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
 * @file tests/gyro-motion-f0/test-gyro-motion-f0.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement tests for gyro-motion-f0 example
 */

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <mpi.h>

#include "config.h"
#include "output-parameters.h"
#include "sapphirepp-logstream.h"
#include "vfp-parameters.h"
#include "vfp-solver.h"


namespace sapinternal
{
  namespace AnalyticSolutionImplementation
  {
    using namespace dealii;
    using namespace sapphirepp;

    template <unsigned int dim>
    class AnalyticSolution : public dealii::Function<dim>
    {
    public:
      AnalyticSolution(const PhysicalParameters &physical_parameters,
                       const unsigned int        exp_order,
                       const double             &mass)
        : Function<dim>((exp_order + 1) * (exp_order + 1))
        , prm{physical_parameters}
        , lms_indices{VFP::PDESystem::create_lms_indices(exp_order)}
        , particle_velocity(mass)
        , scattering_frequency(prm)
      {}



      void
      vector_value(const Point<dim> &point, Vector<double> &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        const double u  = prm.u_sh;
        const double r  = prm.compression_ratio;
        const double p0 = prm.p_inj;

        const double log_p = point[1];
        const double x     = point[0];

        std::vector<double>     values(1);
        std::vector<Point<dim>> points = {point};
        particle_velocity.value_list(points, values);
        const double v = values[0];
        scattering_frequency.value_list(points, values);
        const double nu = values[0];

        f = 0;

        f[0] = std::sqrt(4 * M_PI) * prm.Q / (p0 * u) * 3. * r / (r - 1.) *
               std::exp(-3 * r / (r - 1.) * (log_p - std::log(p0)));

        if (x < 0.)
          {
            f[0] *= std::exp(3 * u * nu / (v * v) * x);
            f[2] = f[0] * (-3. * u / v) / std::sqrt(3);
          }
      }



    private:
      const PhysicalParameters                       prm;
      const std::vector<std::array<unsigned int, 3>> lms_indices;
      const VFP::ParticleVelocity<dim, true>         particle_velocity;
      const VFP::ScatteringFrequency<dim>            scattering_frequency;
    };



    template <unsigned int dim>
    class WeightFunction : public Function<dim>
    {
    public:
      WeightFunction(const PhysicalParameters &physical_parameters,
                     const unsigned int        exp_order)
        : Function<dim>((exp_order + 1) * (exp_order + 1))
        , prm{physical_parameters}
        , lms_indices{VFP::PDESystem::create_lms_indices(exp_order)}
      {}



      void
      vector_value(const Point<dim> &point,
                   Vector<double>   &weight) const override
      {
        AssertDimension(weight.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        const double p0    = prm.p_inj;
        const double log_p = point[1];

        weight = 0;

        if (std::exp(log_p) > p0)
          {
            weight[0] = std::exp(4 * log_p);
            weight[2] = std::exp(4 * log_p);
          }
      }



    private:
      const PhysicalParameters                       prm;
      const std::vector<std::array<unsigned int, 3>> lms_indices;
    };
  } // namespace AnalyticSolutionImplementation
} // namespace sapinternal



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

      saplog.init(2);

      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      double max_L2_error = 1e-10;
      if (argc > 2)
        max_L2_error = std::stod(argv[2]);

      dealii::Timer            timer;
      ParameterHandler         prm;
      VFPParameters<dimension> vfp_parameters;
      PhysicalParameters       physical_parameters;
      Utils::OutputParameters  output_parameters;

      vfp_parameters.declare_parameters(prm);
      physical_parameters.declare_parameters(prm);
      output_parameters.declare_parameters(prm);

      prm.parse_input(parameter_filename);

      vfp_parameters.parse_parameters(prm);
      physical_parameters.parse_parameters(prm);
      output_parameters.parse_parameters(prm);

      timer.start();
      VFPSolver<dimension> vfp_solver(vfp_parameters,
                                      physical_parameters,
                                      output_parameters);
      vfp_solver.run();
      timer.stop();


      saplog << "Calculate analytic solution" << std::endl;
      using namespace sapinternal::AnalyticSolutionImplementation;

      AnalyticSolution<dimension> analytic_solution(
        physical_parameters,
        vfp_parameters.expansion_order,
        vfp_parameters.mass);

      WeightFunction<dimension> weight(physical_parameters,
                                       vfp_parameters.expansion_order);

      const double L2_error =
        vfp_solver.compute_global_error(analytic_solution,
                                        dealii::VectorTools::L2_norm,
                                        dealii::VectorTools::L2_norm,
                                        &weight);
      const double L2_norm =
        vfp_solver.compute_weighted_norm(dealii::VectorTools::L2_norm,
                                         dealii::VectorTools::L2_norm,
                                         &weight);

      PETScWrappers::MPI::Vector analytic_solution_vector;
      analytic_solution_vector.reinit(
        vfp_solver.get_dof_handler().locally_owned_dofs(), MPI_COMM_WORLD);
      dealii::VectorTools::interpolate(vfp_solver.get_dof_handler(),
                                       analytic_solution,
                                       analytic_solution_vector);
      PETScWrappers::MPI::Vector weight_vector;
      weight_vector.reinit(vfp_solver.get_dof_handler().locally_owned_dofs(),
                           MPI_COMM_WORLD);
      dealii::VectorTools::interpolate(vfp_solver.get_dof_handler(),
                                       weight,
                                       weight_vector);

      const unsigned int num_exp_coefficients =
        (vfp_parameters.expansion_order + 1) *
        (vfp_parameters.expansion_order + 1);
      std::vector<std::string> component_names(num_exp_coefficients);
      std::vector<std::string> component_names_weight(num_exp_coefficients);
      const std::vector<std::array<unsigned int, 3>> lms_indices =
        PDESystem::create_lms_indices(num_exp_coefficients);
      for (unsigned int i = 0; i < num_exp_coefficients; ++i)
        {
          const std::array<unsigned int, 3> &lms = lms_indices[i];
          component_names[i] = "analytic_f_" + std::to_string(lms[0]) +
                               std::to_string(lms[1]) + std::to_string(lms[2]);
          component_names_weight[i] = "weight_" + std::to_string(lms[0]) +
                                      std::to_string(lms[1]) +
                                      std::to_string(lms[2]);
        }

      dealii::DataOut<dimension> data_out;
      data_out.attach_dof_handler(vfp_solver.get_dof_handler());
      data_out.add_data_vector(analytic_solution_vector, component_names);
      data_out.add_data_vector(weight_vector, component_names_weight);
      data_out.build_patches(vfp_parameters.polynomial_degree);
      output_parameters.base_file_name = "analytic_solution";
      output_parameters.write_results<dimension>(data_out, 0);

      saplog << "L2_error = " << L2_error << ", L2_norm = " << L2_norm
             << ", rel error = " << L2_error / L2_norm
             << ", CPU/wall time = " << timer.cpu_time() << "/"
             << timer.wall_time() << " s" << std::endl;

      AssertThrow((L2_error / L2_norm) < max_L2_error,
                  dealii::ExcMessage(
                    "L2 error is too large! (" +
                    dealii::Utilities::to_string(L2_error / L2_norm) + " > " +
                    dealii::Utilities::to_string(max_L2_error) + ")"));
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
