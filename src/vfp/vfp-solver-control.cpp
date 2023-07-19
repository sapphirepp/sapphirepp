#include "vfp-solver-control.h"

#include <deal.II/base/patterns.h>

#include <cctype>
#include <sstream>

#include "parameter-flags.h"
#include "parameter-parser.h"

Sapphire::VFPSolverControl::VFPSolverControl(
  const Sapphire::Utils::ParameterParser &prm)
  : expansion_order{prm.vfp_expansion_order}
  , polynomial_degree{prm.vfp_polynomial_degree}
  , theta{prm.vfp_theta}
  , time_step{prm.vfp_time_step}
  , final_time{prm.vfp_final_time}
{
  // Mesh
  // Two diagonally opposite corner points of the grid
  std::stringstream p1_string(prm.vfp_p1);
  std::stringstream p2_string(prm.vfp_p2);
  for (auto [coordinate, i] =
         std::tuple<std::string, unsigned int>{std::string(), 0};
       std::getline(p1_string, coordinate, ',');
       ++i)
    {
      if (i < dim)
        p1[i] = std::stod(coordinate);
    }
  for (auto [coordinate, i] =
         std::tuple<std::string, unsigned int>{std::string(), 0};
       std::getline(p2_string, coordinate, ',');
       ++i)
    {
      if (i < dim)
        p2[i] = std::stod(coordinate);
    }

  // Number of cells
  std::stringstream n_cells_string(prm.vfp_n_cells);
  for (std::string n; std::getline(n_cells_string, n, ',');)
    n_cells.push_back(std::stoi(n));
  n_cells.resize(dim); // shrink to dim of the phase-space

  // Periodicity
  std::string periodicity_string = prm.vfp_periodicity;
  // Remove whitespace
  periodicity_string.erase(std::remove_if(periodicity_string.begin(),
                                          periodicity_string.end(),
                                          [](unsigned char x) {
                                            return std::isspace(x);
                                          }),
                           periodicity_string.end());

  std::stringstream string_stream(periodicity_string);
  for (std::string value; std::getline(string_stream, value, ',');)
    periodicity.push_back((value == "true" ? true : false));
  periodicity.resize(dim);

  // Time stepping method
  if (prm.vfp_time_stepping_method == "Forward Euler")
    {
      theta                = 0.;
      time_stepping_method = Sapphire::Utils::TimeSteppingMethod::forward_euler;
    }
  else if (prm.vfp_time_stepping_method == "Backward Euler")
    {
      theta = 1.;
      time_stepping_method =
        Sapphire::Utils::TimeSteppingMethod::backward_euler;
    }
  else if (prm.vfp_time_stepping_method == "Crank-Nicolson")
    {
      theta = 1. / 2.;
      time_stepping_method =
        Sapphire::Utils::TimeSteppingMethod::crank_nicolson;
    }
  else if (prm.vfp_time_stepping_method == "ERK4")
    time_stepping_method = Sapphire::Utils::TimeSteppingMethod::erk4;
  else if (prm.vfp_time_stepping_method == "LSERK")
    time_stepping_method = Sapphire::Utils::TimeSteppingMethod::lserk;
}

// void
// Sapphire::VFPSolverControl::print_settings(std::ostream &os) const
// {
//   os << "Compile time parameters: " << std::endl;
//   os << "	Dimension Configuration Space: "
//      << VFPSolverControl::dim_configuration_space << "\n";
//   os << "	" << terms << std::endl;
//   os << "Runtime parameters:" << std::endl;
//   parameter_handler.print_parameters(os,
//   dealii::ParameterHandler::ShortText); os << std::endl;
// }

// void
// Sapphire::VFPSolverControl::declare_parameters()
// {
//   parameter_handler.enter_subsection("Mesh");
//   { // NOTE: This is a very strange syntax
//     parameter_handler.declare_entry(
//       "Point 1",
//       "0., 0., 0.",
//       dealii::Patterns::Anything(),
//       "Two diagonally opposite corner points, Point 1 and  Point 2");
//     parameter_handler.declare_entry(
//       "Point 2",
//       "1., 1., 1.",
//       dealii::Patterns::Anything(),
//       "Two diagonally opposite corner points Point 1 and  Point 2");
//     parameter_handler.declare_entry(
//       "Number of cells",
//       "4, 4, 4",
//       dealii::Patterns::Anything(),
//       "Number of cells in each coordinate direction");
//     parameter_handler.declare_entry(
//       "Periodicity",
//       "false, false, false",
//       dealii::Patterns::Anything(),
//       "Periodic boundaries in the three coordinate directions.");
//   }
//   parameter_handler.leave_subsection();

//   parameter_handler.enter_subsection("Time stepping");
//   {
//     parameter_handler.declare_entry(
//       "Method",
//       "Crank-Nicolson",
//       dealii::Patterns::Selection(
//         "Forward Euler|Backward Euler|Crank-Nicolson|ERK4|LSERK"),
//       "The time stepping method.");
//     parameter_handler.declare_entry("Time step size",
//                                     "7.8125e-3",
//                                     dealii::Patterns::Double(),
//                                     "Duration of the simulation.");
//     parameter_handler.declare_entry("Final time",
//                                     "0.4",
//                                     dealii::Patterns::Double(),
//                                     "Duration of the simulation.");
//   }
//   parameter_handler.leave_subsection();

//   parameter_handler.enter_subsection("Expansion");
//   {
//     parameter_handler.declare_entry(
//       "Expansion order",
//       "0",
//       dealii::Patterns::Integer(0),
//       "The order of the expansion of the particle distribution function.");
//   }
//   parameter_handler.leave_subsection();

//   parameter_handler.enter_subsection("Finite element");
//   {
//     parameter_handler.declare_entry("Polynomial degree",
//                                     "1",
//                                     dealii::Patterns::Integer(0),
//                                     "The degree of the shape functions (i.e.
//                                     " "the polynomials) of the finite
//                                     element.");
//   }
//   parameter_handler.leave_subsection();
//   parameter_handler.enter_subsection("Output");
//   {
//     parameter_handler.declare_entry(
//       "Results folder",
//       "./results",
//       dealii::Patterns::Anything(),
//       "Path to the folder in which the simulation results will be stored. "
//       "Without a trailing slash.");
//     parameter_handler.declare_entry(
//       "Simulation identifier",
//       "001",
//       dealii::Patterns::Anything(),
//       "Name of the simulation run. It will be used to create a subdirectory "
//       "in the results folder.");
//     parameter_handler.declare_entry(
//       "Format",
//       "vtu",
//       dealii::Patterns::Selection("vtu|hdf5"),
//       "The format in which the simulation output will be stored.");
//     parameter_handler.declare_entry("Output frequency",
//                                     "1",
//                                     dealii::Patterns::Integer(0),
//                                     "The frequence at which output files will
//                                     " "be written. (In units of time
//                                     steps)");
//   }
//   parameter_handler.leave_subsection();
// }

// void
// Sapphire::VFPSolverControl::parse_parameters()
// {
//   parameter_handler.parse_input(parameter_file);
// }

// void
// Sapphire::VFPSolverControl::get_parameters()
// {
//   parameter_handler.enter_subsection("Mesh");
//   {
//     // Two diagonally opposite corner points of the grid
//     std::stringstream p1_string(parameter_handler.get("Point 1"));
//     std::stringstream p2_string(parameter_handler.get("Point 2"));
//     for (auto [coordinate, i] =
//            std::tuple<std::string, unsigned int>{std::string(), 0};
//          std::getline(p1_string, coordinate, ',');
//          ++i)
//       {
//         if (i < dim)
//           p1[i] = std::stod(coordinate);
//       }
//     for (auto [coordinate, i] =
//            std::tuple<std::string, unsigned int>{std::string(), 0};
//          std::getline(p2_string, coordinate, ',');
//          ++i)
//       {
//         if (i < dim)
//           p2[i] = std::stod(coordinate);
//       }

//     // Number of cells
//     std::stringstream n_cells_string(parameter_handler.get("Number of
//     cells")); for (std::string n; std::getline(n_cells_string, n, ',');)
//       n_cells.push_back(std::stoi(n));
//     n_cells.resize(dim); // shrink to dim of the phase-space

//     // Periodicity
//     std::string periodicity_string = parameter_handler.get("Periodicity");
//     // Remove whitespace
//     periodicity_string.erase(std::remove_if(periodicity_string.begin(),
//                                             periodicity_string.end(),
//                                             [](unsigned char x) {
//                                               return std::isspace(x);
//                                             }),
//                              periodicity_string.end());

//     std::stringstream string_stream(periodicity_string);
//     for (std::string value; std::getline(string_stream, value, ',');)
//       periodicity.push_back((value == "true" ? true : false));
//     periodicity.resize(dim);
//   }
//   parameter_handler.leave_subsection();

//   parameter_handler.enter_subsection("Time stepping");
//   {
//     time_stepping_method = parameter_handler.get("Method");
//     if (time_stepping_method == "Forward Euler")
//       theta = 0.;
//     else if (time_stepping_method == "Backward Euler")
//       theta = 1.;
//     else if (time_stepping_method == "Crank-Nicolson")
//       theta = 1. / 2.;
//     else
//       theta = 1. / 2.; // default
//     time_step  = parameter_handler.get_double("Time step size");
//     final_time = parameter_handler.get_double("Final time");
//   }
//   parameter_handler.leave_subsection();

//   parameter_handler.enter_subsection("Expansion");
//   {
//     expansion_order = parameter_handler.get_integer("Expansion order");
//   }
//   parameter_handler.leave_subsection();

//   parameter_handler.enter_subsection("Finite element");
//   {
//     polynomial_degree = parameter_handler.get_integer("Polynomial degree");
//   }
//   parameter_handler.leave_subsection();
//   parameter_handler.enter_subsection("Output");
//   {
//     results_path     = parameter_handler.get("Results folder");
//     simulation_id    = parameter_handler.get("Simulation identifier");
//     format           = parameter_handler.get("Format");
//     output_frequency = parameter_handler.get_integer("Output frequency");
//   }
// }
