#include <deal.II/base/mpi.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <mpi.h>

#include "output-module.h"
#include "vfp-equation-solver.h"


/**
 * @brief Creates a 2d mesh around a shock
 * @param filename File where the mesh is saved
 * @param physical_properties Physical properties with defined grid parameters
 *
 * Create a 2d mesh with a shock at x = 0. In the shock region the mesh is
 * equidistant with a constant step size. Outside the shock region the step size
 * increases with a scaling of the step size:
 * \f[ \Delta x_i = (scaling delta_x)^i \cdot \Delta x_{\rm shock} \f]
 * The grid in momentum direction is equidistant.
 */
void
make_grid_shock(const std::string                  &filename,
                const Sapphire::PhysicalProperties &physical_properties)
{
  const unsigned int n_cells_downstream =
    physical_properties.n_cells_downstream;
  const unsigned int n_cells_upstream = physical_properties.n_cells_upstream;
  const unsigned int n_cells_shock    = physical_properties.n_cells_shock;
  const double       shock_width      = physical_properties.shock_width;
  const double       scaling_delta_x  = physical_properties.scaling_delta_x;
  const double       log_p_min        = std::log(physical_properties.p_min);
  const double       log_p_max        = std::log(physical_properties.p_max);
  const unsigned int n_cells_p        = physical_properties.n_cells_p;

  using namespace dealii;
  using namespace Sapphire;
  LogStream::Prefix p("ShockGirdGenerator", Sapphire::saplog);
  saplog << "Generate grid \t[" << Utilities::System::get_time() << "]"
         << std::endl;
  parallel::shared::Triangulation<2> triangulation(MPI_COMM_WORLD);

  std::vector<std::vector<double>> step_sizes{
    std::vector<double>(n_cells_downstream + n_cells_upstream +
                        2 * n_cells_shock),
    std::vector<double>(n_cells_p)};

  // equidistant momentum mesh
  const double delta_h_p = (log_p_max - log_p_min) / n_cells_p;
  std::fill(step_sizes[1].begin(), step_sizes[1].end(), delta_h_p);

  // configuration space x
  // TODO: Decide if we want to use sinh or scaling_delta_x
  const double step_size_shock = shock_width / (2 * n_cells_shock);
  //   2 * shock_width / n_cells_shock; // Old wrong calculation
  // const double start_sinh = std::asinh(step_size_shock);

  double length_upstream = 0.;
  for (unsigned int i = 0; i < n_cells_upstream; ++i)
    {
      step_sizes[0][n_cells_upstream - 1 - i] =
        std::pow(scaling_delta_x, i) * step_size_shock;
      // step_sizes[0][n_cells_upstream - 1 - i] =
      //   std::sinh(i * scaling_delta_x + start_sinh);
      length_upstream += step_sizes[0][n_cells_upstream - 1 - i];
    }

  for (unsigned int i = 0; i < 2 * n_cells_shock; ++i)
    step_sizes[0][n_cells_upstream + i] = step_size_shock;

  double length_downstream = 0.;
  for (unsigned int i = 0; i < n_cells_downstream; ++i)
    {
      step_sizes[0][n_cells_upstream + 2 * n_cells_shock + i] =
        std::pow(scaling_delta_x, i) * step_size_shock;
      // step_sizes[0][n_cells_upstream + 2 * n_cells_shock + i] =
      //   std::sinh(i * scaling_delta_x + start_sinh);
      length_downstream +=
        step_sizes[0][n_cells_upstream + 2 * n_cells_shock + i];
    }

  length_upstream += step_size_shock * n_cells_shock;
  length_downstream += step_size_shock * n_cells_shock;

  Point<2> p1{-length_upstream, log_p_min};
  Point<2> p2{+length_downstream, log_p_max};

  GridGenerator::subdivided_hyper_rectangle(
    triangulation, step_sizes, p1, p2, true);

  LogStream::Prefix prefix2("GridInfo", saplog);
  saplog << "length_upstream=" << length_upstream
         << ", length_downstream=" << length_downstream
         << ", step_size_shock=" << step_size_shock << ", #cells in x="
         << n_cells_downstream + n_cells_upstream + 2 * n_cells_shock
         << ", #active cells=" << triangulation.n_global_active_cells()
         << std::endl;

  // x_max is given by the scale width kappa, use x_max = 5* kappa
  // we want to achieve upstream > 5 * 1/kappa
  // downstream can be shorter
  const double alpha = physical_properties.alpha;
  const double B_0   = physical_properties.B_0;
  const double u_sh  = physical_properties.u_sh;
  const double kappa = 3. * alpha * B_0 / std::exp(log_p_max) * u_sh;
  saplog << "1/kappa=" << 1. / kappa
         << ", downstream*kappa=" << length_downstream * kappa
         << ", upstream*kappa=" << length_upstream * kappa << std::endl;
  // Assert(length_upstream * kappa > 5.0,
  //        ExcMessage("length upstream too short"));
  // Assert(length_downstream * kappa > 2.0,
  //        ExcMessage("length downstream too short"));

  saplog << "Write grid to \"" << filename << "\"" << std::endl;

  GridOutFlags::Ucd ucd_flags;
  ucd_flags.write_faces = true;
  ucd_flags.write_lines = true;

  GridOut grid_out;
  grid_out.set_flags(ucd_flags);
  std::ofstream output(filename);
  grid_out.write_ucd(triangulation, output);
}


int
main(int argc, char *argv[])
{
  try
    {
      using namespace Sapphire;
      using namespace VFP;
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      saplog.depth_console(2);
      const unsigned int mpi_rank =
        Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      if (mpi_rank > 0)
        {
          saplog.depth_console(0);
          char buffer[4];
          std::snprintf(buffer, 4, "%03d", mpi_rank);
          saplog.push("mpi" + std::string(buffer));
        }

      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      int mpi_size = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      saplog << "Start oblique_shock with parameter file \""
             << parameter_filename << "\" on " << mpi_size << " processor(s) ["
             << Utilities::System::get_date() << " "
             << Utilities::System::get_time() << "]" << std::endl;

      saplog.push("main");
      ParameterHandler       prm;
      VFPSolverControl       vfp_solver_control;
      PhysicalProperties     physical_properties;
      Utils::OutputModule<2> output_module;

      vfp_solver_control.declare_parameters(prm);
      physical_properties.declare_parameters(prm);
      output_module.declare_parameters(prm);

      prm.parse_input(parameter_filename);

      vfp_solver_control.parse_parameters(prm);
      physical_properties.parse_parameters(prm);
      output_module.parse_parameters(prm);

      output_module.init(prm);

      saplog.depth_console(4);
      make_grid_shock(vfp_solver_control.grid_file, physical_properties);
      saplog.depth_console(2);

      saplog.pop();
      VFPEquationSolver vfp_equation_solver(vfp_solver_control,
                                            physical_properties,
                                            output_module);
      vfp_equation_solver.run();

      // TODO: Calculate error
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
