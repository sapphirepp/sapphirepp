#include "vfp-solver-control.h"

#include <deal.II/base/patterns.h>

#include <cctype>
#include <sstream>

#include "parameter-flags.h"
#include "parameter-parser.h"

Sapphire::VFP::VFPSolverControl::VFPSolverControl(
  const Sapphire::Utils::ParameterParser &prm)
  : expansion_order{prm.vfp_expansion_order}
  , polynomial_degree{prm.vfp_polynomial_degree}
  , theta{prm.vfp_theta}
  , time_step{prm.vfp_time_step}
  , grid_type(prm.vfp_grid_type)
  , grid_file(prm.vfp_grid_file)
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
