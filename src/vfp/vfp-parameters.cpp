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
 * @file vfp-parameters.cpp
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement @ref sapphirepp::VFP::VFPParameters
 */

#include "vfp-parameters.h"

#include <deal.II/base/patterns.h>

#include <filesystem>
#include <sstream>

#include "sapphirepp-logstream.h"



template <unsigned int dim>
sapphirepp::VFP::VFPParameters<dim>::VFPParameters() = default;



template <unsigned int dim>
void
sapphirepp::VFP::VFPParameters<dim>::declare_parameters(ParameterHandler &prm)
{
  LogStream::Prefix pre1("Startup", saplog);
  LogStream::Prefix pre2("VFPParameters", saplog);
  saplog << "Declaring parameters" << std::endl;
  prm.enter_subsection("VFP");


  prm.enter_subsection("Mesh");
  {
    prm.declare_entry("Grid type",
                      "Shock grid",
                      Patterns::Selection("Hypercube|Shock grid|File"),
                      "The type of the grid. Can either be created by the "
                      "program or read from a file");

    prm.declare_entry("Point 1",
                      "-2., -2., -2.",
                      Patterns::Anything(),
                      "Two diagonally opposite corner points, "
                      "Point 1 and  Point 2");
    prm.declare_entry("Point 2",
                      "2., 2., 2.",
                      Patterns::Anything(),
                      "Two diagonally opposite corner points, "
                      "Point 1 and  Point 2");
    prm.declare_entry("Number of cells",
                      "32, 32, 32",
                      Patterns::Anything(),
                      "Number of cells in each coordinate direction");

    prm.declare_entry("Number of shock cells",
                      "8",
                      Patterns::Integer(1),
                      "Number of cells along x-direction in the shock (for "
                      "Shock grid)");
    prm.declare_entry("Shock width",
                      "0.01",
                      Patterns::Double(0),
                      "Width of the shock in x-direction (for Shock grid)");
    prm.declare_entry("Scaling factor shock",
                      "1.5",
                      Patterns::Double(1),
                      "Scaling factor for the cell width in x-direction "
                      "(for Shock grid)");

    prm.declare_entry("File name",
                      "",
                      Patterns::Anything(),
                      "The file containing the grid (only for "
                      "Grid type = File)");

    prm.enter_subsection("Boundary conditions");
    {
      const auto boundary_pattern =
        Patterns::Selection("continuous|zero inflow|periodic");
      prm.declare_entry("lower x",
                        "continuous",
                        boundary_pattern,
                        "Boundary condition at the lower x boundary.");
      prm.declare_entry("upper x",
                        "continuous",
                        boundary_pattern,
                        "Boundary condition at the upper x boundary.");
      prm.declare_entry("lower y",
                        "continuous",
                        boundary_pattern,
                        "Boundary condition at the lower y boundary.");
      prm.declare_entry("upper y",
                        "continuous",
                        boundary_pattern,
                        "Boundary condition at the upper y boundary.");
      prm.declare_entry("lower z",
                        "continuous",
                        boundary_pattern,
                        "Boundary condition at the lower z boundary.");
      prm.declare_entry("upper z",
                        "continuous",
                        boundary_pattern,
                        "Boundary condition at the upper z boundary.");
      prm.declare_entry("lower p",
                        "continuous",
                        boundary_pattern,
                        "Boundary condition at the lower p boundary.");
      prm.declare_entry("upper p",
                        "continuous",
                        boundary_pattern,
                        "Boundary condition at the upper p boundary.");
    }
    prm.leave_subsection();
  } // Mesh
  prm.leave_subsection();


  prm.enter_subsection("Time stepping");
  {
    prm.declare_entry("Method",
                      "CN",
                      /** @todo "LSERK" is not working, so we exclude it */
                      // Patterns::Selection("FE|BE|CN|ERK4|LSERK"),
                      Patterns::Selection("FE|BE|CN|ERK4"),
                      "The time stepping method.");
    prm.declare_entry("Time step size",
                      "1.0",
                      Patterns::Double(0),
                      "Time step for the simulation in dimensionless units.");
    prm.declare_entry("Final time",
                      "200",
                      Patterns::Double(0),
                      "End time for the simulation in dimensionless units.");
  } // Time stepping
  prm.leave_subsection();


  prm.enter_subsection("Expansion");
  {
    prm.declare_entry("Expansion order",
                      "1",
                      Patterns::Integer(0),
                      "The maximum order of the expansion of the particle "
                      "distribution function.");
  } // Expansion
  prm.leave_subsection();


  prm.enter_subsection("Finite element");
  {
    prm.declare_entry("Polynomial degree",
                      "1",
                      Patterns::Integer(0),
                      "The degree of the shape functions (i.e. the "
                      "polynomials) of the finite element.");
  } // Finite element
  prm.leave_subsection();


  prm.enter_subsection("Particle properties");
  {
    prm.declare_entry("Mass",
                      "1.",
                      Patterns::Double(0),
                      "Mass of the particles in dimensionless units.");
    prm.declare_entry("Charge",
                      "1.",
                      Patterns::Double(),
                      "The charge of the particles in dimensionless units.");
  } // Particle properties
  prm.leave_subsection();


  prm.enter_subsection("TransportOnly");
  {
    prm.declare_entry("Gamma",
                      "3.",
                      Patterns::Double(1.),
                      "The Lorentz factor of the particles. "
                      "Only has to be specified in the transport-only (i.e. "
                      "p-independent) case.");
  } // TransportOnly
  prm.leave_subsection();


  prm.leave_subsection();
}



template <unsigned int dim>
void
sapphirepp::VFP::VFPParameters<dim>::parse_parameters(ParameterHandler &prm)
{
  LogStream::Prefix pre1("Startup", saplog);
  LogStream::Prefix pre2("VFPParameters", saplog);
  saplog << "Parsing parameters" << std::endl;
  std::string s;
  prm.enter_subsection("VFP");


  prm.enter_subsection("Mesh");
  {
    s = prm.get("Grid type");
    if (s == "Hypercube")
      grid_type = GridType::hypercube;
    else if (s == "Shock grid")
      grid_type = GridType::shock;
    else if (s == "File")
      grid_type = GridType::file;
    else
      Assert(false, ExcNotImplemented());

    // Two diagonally opposite corner points of the grid
    unsigned int i = 0;
    s              = prm.get("Point 1");
    std::stringstream p1_string(s);
    for (std::string coordinate; std::getline(p1_string, coordinate, ','); ++i)
      {
        if (i < dim)
          p1[i] = std::stod(coordinate);
      }
    AssertThrow(dim <= i,
                ExcMessage(
                  "Point 1 does not specify coordinate in each dimension. "
                  "Please enter the coordinates of the lower left corner "
                  "of the grid: \n"
                  "\tset Point 1 = x1 (, y1) (, z1) (, p1) \n"
                  "You entered: " +
                  s));
    if (i != dim)
      saplog << "WARNING: Point 1 specification does not match dimension!"
             << std::endl;

    i = 0;
    s = prm.get("Point 2");
    std::stringstream p2_string(s);
    for (std::string coordinate; std::getline(p2_string, coordinate, ','); ++i)
      {
        if (i < dim)
          p2[i] = std::stod(coordinate);
      }
    AssertThrow(dim <= i,
                ExcMessage(
                  "Point 2 does not specify coordinate in each dimension. "
                  "Please enter the coordinates of the lower left corner "
                  "of the grid: \n"
                  "\tset Point 2 = x1 (, y1) (, z1) (, p1) \n"
                  "You entered: " +
                  s));
    if (i != dim)
      saplog << "WARNING: Point 2 specification does not match dimension!"
             << std::endl;

    // Number of cells
    s = prm.get("Number of cells");
    std::stringstream n_cells_string(s);
    for (std::string n; std::getline(n_cells_string, n, ',');)
      n_cells.push_back(static_cast<unsigned int>(std::stoi(n)));
    AssertThrow(dim <= n_cells.size(),
                ExcMessage(
                  "Number of cells does not specify value in each dimension."
                  "Please enter the number of cells in each coordinate: \n"
                  "\tset Number of cells = Nx (, Ny) (, Nz) (, Np) \n"
                  "You entered: " +
                  s));
    if (n_cells.size() != dim)
      saplog
        << "WARNING: Number of cells specification does not match dimension!"
        << std::endl;
    n_cells.resize(dim);

    n_cells_shock =
      static_cast<unsigned int>(prm.get_integer("Number of shock cells"));
    shock_width          = prm.get_double("Shock width");
    scaling_factor_shock = prm.get_double("Scaling factor shock");

    grid_file = prm.get("File name");

    prm.enter_subsection("Boundary conditions");
    {
      boundary_conditions.resize(2 * dim);

      for (unsigned int boundary_id = 0; boundary_id < 2 * dim; ++boundary_id)
        {
          std::string entry = "";
          if (boundary_id % 2 == 0)
            entry = "lower ";
          else
            entry = "upper ";

          if (boundary_id / 2 == dim_cs)
            entry += "p";
          else if (boundary_id / 2 == 0)
            entry += "x";
          else if (boundary_id / 2 == 1)
            entry += "y";
          else if (boundary_id / 2 == 2)
            entry += "z";
          else
            Assert(false, ExcNotImplemented());

          s = prm.get(entry);
          if (s == "continuous")
            boundary_conditions[boundary_id] =
              VFP::BoundaryConditions::continuous;
          else if (s == "zero inflow")
            boundary_conditions[boundary_id] =
              VFP::BoundaryConditions::zero_inflow;
          else if (s == "periodic")
            boundary_conditions[boundary_id] =
              VFP::BoundaryConditions::periodic;
          else
            Assert(false, ExcNotImplemented());
        }
    }
    prm.leave_subsection();
  } // Mesh
  prm.leave_subsection();


  prm.enter_subsection("Time stepping");
  {
    s = prm.get("Method");
    if (s == "FE")
      {
        theta                = 0.;
        time_stepping_method = TimeSteppingMethod::forward_euler;
      }
    else if (s == "BE")
      {
        theta                = 1.;
        time_stepping_method = TimeSteppingMethod::backward_euler;
      }
    else if (s == "CN")
      {
        theta                = 0.5;
        time_stepping_method = TimeSteppingMethod::crank_nicolson;
      }
    else if (s == "ERK4")
      time_stepping_method = TimeSteppingMethod::erk4;
    else if (s == "LSERK")
      time_stepping_method = TimeSteppingMethod::lserk;
    else
      Assert(false, ExcNotImplemented());

    time_step  = prm.get_double("Time step size");
    final_time = prm.get_double("Final time");
  } // Time stepping
  prm.leave_subsection();


  prm.enter_subsection("Expansion");
  {
    expansion_order =
      static_cast<unsigned int>(prm.get_integer("Expansion order"));
  } // Expansion
  prm.leave_subsection();


  prm.enter_subsection("Finite element");
  {
    polynomial_degree =
      static_cast<unsigned int>(prm.get_integer("Polynomial degree"));
  } // Finite element
  prm.leave_subsection();


  prm.enter_subsection("Particle properties");
  {
    mass   = prm.get_double("Mass");
    charge = prm.get_double("Charge");
  } // Particle properties
  prm.leave_subsection();


  prm.enter_subsection("TransportOnly");
  {
    gamma    = prm.get_double("Gamma");
    velocity = std::sqrt(1 - 1 / (gamma * gamma));
  } // TransportOnly
  prm.leave_subsection();


  prm.leave_subsection();
}



// Explicit instantiation
template class sapphirepp::VFP::VFPParameters<1>;
template class sapphirepp::VFP::VFPParameters<2>;
template class sapphirepp::VFP::VFPParameters<3>;
