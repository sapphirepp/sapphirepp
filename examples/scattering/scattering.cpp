/////////////////////////////////////////////////
/// @file examples/scattering/scattering.cpp
/// @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
/// @brief Implements the main for the scattering example
/// @date 2023-10-30
///
/// @copyright Copyright (c) 2023
///
/////////////////////////////////////////////////

/// @page scattering
/// @subsection main scattering.cpp
/// @dontinclude scattering/scattering.cpp
/// In this file, we will define the `main` function for the scattering example.
/// First, we include the header files we need:
#include <deal.II/base/mpi.h>

#include <mpi.h>

#include "config.h"
#include "output-module.h"
#include "vfp-equation-solver.h"

/// @page scattering
/// @skip include
/// @until CODEBLOCK_DELIMITER
/// Then we define the `main` function:
int
main(int argc, char *argv[])
{
  /// @page scattering
  /// @until CODEBLOCK_DELIMITER
  /// We will run the program inside a try-catch block to catch any exceptions
  /// and output them.
  try
    {
      /// @page scattering
      /// @until CODEBLOCK_DELIMITER
      /// First, we need to initialize MPI for parallel runs:
      using namespace Sapphire;
      using namespace VFP;
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);

      /// @page scattering
      /// @until CODEBLOCK_DELIMITER
      /// We will use the `Sapphire::Log` class to output information to the
      /// console. The `depth_console()` functions allows to specify how
      /// detailed the output should be. We will set it to 2, which results in
      /// one line of output per time-step. For parallel runs, we will turn of
      /// the output of all but the first process.
      saplog.depth_console(2);
      saplog.depth_console(-1);
      const unsigned int mpi_rank =
        dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      if (mpi_rank > 0)
        {
          saplog.depth_console(0);
          saplog.push("mpi" + dealii::Utilities::to_string(mpi_rank, 3));
        }

      /// @page scattering
      /// @until CODEBLOCK_DELIMITER
      /// The program allows to specify the parameter file as a command line
      /// argument. In case no argument is given, it will default to the file
      /// `parameter.prm`. When starting the program, we output the
      /// configuration specified by command line arguments.
      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      int mpi_size = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      saplog << "Start example scattering with parameter file \""
             << parameter_filename << "\" on " << mpi_size << " processor(s) ["
             << dealii::Utilities::System::get_date() << " "
             << dealii::Utilities::System::get_time() << "]" << std::endl;

      /// @page scattering
      /// @until CODEBLOCK_DELIMITER
      /// Next, we create all object that declare runtime parameter.
      /// We prefix these lines with `main` in the log stream, to silence the
      /// initialisation process in a run with low verbosity.
      saplog.push("main");
      ParameterHandler               prm;
      VFPSolverControl<dimension>    vfp_solver_control;
      PhysicalProperties             physical_properties;
      Utils::OutputModule<dimension> output_module;

      /// @page scattering
      /// @until CODEBLOCK_DELIMITER
      /// We now declare all parameters that the ParameterHandler should expect.
      /// @note All parameters have to be declared before we can parse the
      /// parameter file.
      vfp_solver_control.declare_parameters(prm);
      physical_properties.declare_parameters(prm);
      output_module.declare_parameters(prm);

      /// @page scattering
      /// @until CODEBLOCK_DELIMITER
      /// After declaring the parameters, we can parse the parameter file,
      prm.parse_input(parameter_filename);

      /// @page scattering
      /// @until CODEBLOCK_DELIMITER
      /// and then parse the parameters to objects.
      vfp_solver_control.parse_parameters(prm);
      physical_properties.parse_parameters(prm);
      output_module.parse_parameters(prm);

      saplog.pop(); // pop "main" prefix

      /// @page scattering
      /// @until CODEBLOCK_DELIMITER
      /// Finally, we can create the VFP equation solver and run it.
      VFPEquationSolver<dimension> vfp_equation_solver(vfp_solver_control,
                                                       physical_properties,
                                                       output_module);
      vfp_equation_solver.run();

      /// @page scattering
      /// @until CODEBLOCK_DELIMITER
      /// After the simulation ended, we can compute the error by comparing to
      /// the analytic solution (implemented in `InitialValueFunction`).
      InitialValueFunction<dimension> analytic_solution(
        physical_properties, vfp_solver_control.expansion_order);
      analytic_solution.set_time(vfp_solver_control.final_time);

      const double L2_error =
        vfp_equation_solver.compute_global_error(analytic_solution,
                                                 VectorTools::L2_norm,
                                                 VectorTools::L2_norm);

      saplog << "L2 error: " << L2_error << std::endl;
    }
  /// @page scattering
  /// @until CODEBLOCK_DELIMITER
  /// In case an exception is thrown, we will catch it, output it ot std::err
  /// and return with an error code.
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

/// @page scattering
/// @until END_OF_FILE
///
/// This concludes the `main` function for the scattering example.
