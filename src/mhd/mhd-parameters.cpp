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
 * @file mhd-parameters.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement @ref sapphirepp::MHD::MHDParameters
 */

#include "mhd-parameters.h"

#include <deal.II/base/patterns.h>

#include <filesystem>
#include <sstream>

#include "sapphirepp-logstream.h"



template <unsigned int dim>
sapphirepp::MHD::MHDParameters<dim>::MHDParameters() = default;



template <unsigned int dim>
void
sapphirepp::MHD::MHDParameters<dim>::declare_parameters(ParameterHandler &prm)
{
  LogStream::Prefix pre1("Startup", saplog);
  LogStream::Prefix pre2("MHDParameters", saplog);
  saplog << "Declaring parameters" << std::endl;
  prm.enter_subsection("MHD");


  prm.enter_subsection("Mesh");
  {
    prm.declare_entry("Grid type",
                      "Hypercube",
                      Patterns::Selection("Hypercube|File"),
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
                      "8, 8, 8",
                      Patterns::Anything(),
                      "Number of cells in each coordinate direction");

    prm.declare_entry("File name",
                      "",
                      Patterns::Anything(),
                      "The file containing the grid (only for "
                      "Grid type = File)");

    prm.enter_subsection("Boundary conditions");
    {
      const auto boundary_pattern = Patterns::Selection("periodic");
      prm.declare_entry("lower x",
                        "periodic",
                        boundary_pattern,
                        "Boundary condition at the lower x boundary.");
      prm.declare_entry("upper x",
                        "periodic",
                        boundary_pattern,
                        "Boundary condition at the upper x boundary.");
      prm.declare_entry("lower y",
                        "periodic",
                        boundary_pattern,
                        "Boundary condition at the lower y boundary.");
      prm.declare_entry("upper y",
                        "periodic",
                        boundary_pattern,
                        "Boundary condition at the upper y boundary.");
      prm.declare_entry("lower z",
                        "periodic",
                        boundary_pattern,
                        "Boundary condition at the lower z boundary.");
      prm.declare_entry("upper z",
                        "periodic",
                        boundary_pattern,
                        "Boundary condition at the upper z boundary.");
    }
    prm.leave_subsection();
  } // Mesh
  prm.leave_subsection();


  prm.enter_subsection("Time stepping");
  {
    prm.declare_entry("Method",
                      "FE",
                      Patterns::Selection("FE|ERK2|ERK4"),
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


  prm.enter_subsection("Finite element");
  {
    prm.declare_entry("Polynomial degree",
                      "1",
                      Patterns::Integer(0),
                      "The degree of the shape functions (i.e. the "
                      "polynomials) of the finite element.");
  } // Finite element
  prm.leave_subsection();


  prm.enter_subsection("Plasma properties");
  {
    prm.declare_entry("Adiabatic index",
                      "1.666666",
                      Patterns::Double(1),
                      "Adiabatic index (gamma) of the plasma.");
  } // Plasma properties
  prm.leave_subsection();


  prm.leave_subsection();
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDParameters<dim>::parse_parameters(ParameterHandler &prm)
{
  LogStream::Prefix pre1("Startup", saplog);
  LogStream::Prefix pre2("MHDParameters", saplog);
  saplog << "Parsing parameters" << std::endl;
  std::string s;
  prm.enter_subsection("MHD");


  prm.enter_subsection("Mesh");
  {
    s = prm.get("Grid type");
    if (s == "Hypercube")
      grid_type = GridTypeMHD::hypercube;
    else if (s == "File")
      grid_type = GridTypeMHD::file;
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
                  "\tset Point 1 = x1 (, y1) (, z1) \n"
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
                  "\tset Point 2 = x1 (, y1) (, z1) \n"
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
                  "\tset Number of cells = Nx (, Ny) (, Nz) \n"
                  "You entered: " +
                  s));
    if (n_cells.size() != dim)
      saplog
        << "WARNING: Number of cells specification does not match dimension!"
        << std::endl;
    n_cells.resize(dim);

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

          if (boundary_id / 2 == 0)
            entry += "x";
          else if (boundary_id / 2 == 1)
            entry += "y";
          else if (boundary_id / 2 == 2)
            entry += "z";
          else
            Assert(false, ExcNotImplemented());

          s = prm.get(entry);
          if (s == "zero inflow")
            boundary_conditions[boundary_id] =
              MHD::BoundaryConditionsMHD::zero_inflow;
          else if (s == "periodic")
            boundary_conditions[boundary_id] =
              MHD::BoundaryConditionsMHD::periodic;
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
      time_stepping_method = TimeSteppingMethodMHD::forward_euler;
    else if (s == "ERK2")
      time_stepping_method = TimeSteppingMethodMHD::erk2;
    else if (s == "ERK4")
      time_stepping_method = TimeSteppingMethodMHD::erk4;
    else
      Assert(false, ExcNotImplemented());

    time_step  = prm.get_double("Time step size");
    final_time = prm.get_double("Final time");
  } // Time stepping
  prm.leave_subsection();


  prm.enter_subsection("Finite element");
  {
    polynomial_degree =
      static_cast<unsigned int>(prm.get_integer("Polynomial degree"));
  } // Finite element
  prm.leave_subsection();


  prm.enter_subsection("Plasma properties");
  {
    adiabatic_index = prm.get_double("Adiabatic index");
  } // Plasma properties
  prm.leave_subsection();


  prm.leave_subsection();
}



// Explicit instantiation
template class sapphirepp::MHD::MHDParameters<1>;
template class sapphirepp::MHD::MHDParameters<2>;
template class sapphirepp::MHD::MHDParameters<3>;
