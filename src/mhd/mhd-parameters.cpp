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
#include <deal.II/base/utilities.h>

#include <filesystem>
#include <limits>

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

  const auto pattern_point = Patterns::List(Patterns::Double(), dim, dim, ",");


  prm.enter_subsection("Mesh");
  {
    prm.declare_entry("Grid type",
                      "Hypercube",
                      Patterns::Selection("Hypercube|File"),
                      "The type of the grid. Can either be created by the "
                      "program or read from a file");

    prm.declare_entry("Point 1",
                      (dim == 1) ? "-2." :
                      (dim == 2) ? "-2., -2." :
                                   "-2., -2., -2.",
                      pattern_point,
                      "Two diagonally opposite corner points, "
                      "Point 1 and  Point 2");
    prm.declare_entry("Point 2",
                      (dim == 1) ? "2." :
                      (dim == 2) ? "2., 2." :
                                   "2., 2., 2.",
                      pattern_point,
                      "Two diagonally opposite corner points, "
                      "Point 1 and  Point 2");
    prm.declare_entry("Number of cells",
                      (dim == 1) ? "32" :
                      (dim == 2) ? "32, 32" :
                                   "32, 32, 32",
                      Patterns::List(Patterns::Integer(1), dim, dim, ","),
                      "Number of cells in each coordinate direction");

    prm.declare_entry("File name",
                      "",
                      Patterns::Anything(),
                      "The file containing the grid (only for "
                      "Grid type = File)");

    prm.enter_subsection("Boundary conditions");
    {
      const auto boundary_pattern = Patterns::Selection("zero inflow|periodic");
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

    prm.declare_entry("Courant number",
                      "0.8",
                      Patterns::Double(0),
                      "Courant number/CFL number "
                      "to infer time step from CFL condition.");
    /** @todo Rename parameter to Maximum time step */
    prm.declare_entry("Time step size",
                      "0.",
                      Patterns::Double(0),
                      "Maximum time step for the simulation "
                      "in dimensionless units. "
                      "Set to zero to infer time step from CFL condition.");
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
                      Patterns::Double(1, 2),
                      "Adiabatic index (gamma) of the plasma.");
  } // Plasma properties
  prm.leave_subsection();


  prm.enter_subsection("Numerical parameters");
  {
    prm.add_parameter(
      "mhd_floor",
      mhd_floor,
      "This section collects numerical parameters for the MHD solver. \n"
      "!DO NOT CHANGE unless you are aware what the parameters do! \n\n"
      "Floor value for pressure, density and energy.",
      Patterns::Double(0.));

    prm.add_parameter("indicator_threshold",
                      indicator_threshold,
                      "Threshold for KXRCF shock indicator.",
                      Patterns::Double(0.));

    prm.add_parameter("minmod_threshold",
                      minmod_threshold,
                      "minmod threshold parameter $M$ for slope limiter.",
                      Patterns::Double(0.));
    prm.add_parameter("minmod_beta",
                      minmod_beta,
                      "minmod limiter parameter $\\beta$ for slope limiter.",
                      Patterns::Double(0.));
  } // Numerical parameters
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
    Patterns::Tools::to_value(prm.get("Point 1"), p1);
    Patterns::Tools::to_value(prm.get("Point 2"), p2);

    // Number of cells
    Patterns::Tools::to_value(prm.get("Number of cells"), n_cells);

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

    courant_number = prm.get_double("Courant number");
    max_time_step  = prm.get_double("Time step size");
    if (max_time_step == 0.)
      max_time_step = std::numeric_limits<double>::max();
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
