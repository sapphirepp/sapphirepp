#include "parameter-parser.h"

#include <iostream>

#include "parameter-flags.h"
#include "sapphire-logstream.h"

Sapphire::Utils::ParameterParser::ParameterParser(
  const std::string &prm_file_name)
  : prm_file_name(prm_file_name)
{
  LogStream::Prefix p("ParameterParser", saplog);
  declare_parameters();
  saplog << "Parsing input file \"" << prm_file_name << "\"" << std::endl;
  prm.parse_input(prm_file_name);
  parse_parameters();
  prm.log_parameters(saplog);
}

void
Sapphire::Utils::ParameterParser::write_parameters(
  const std::string &filename) const
{
  LogStream::Prefix p("ParameterParser", saplog);
  saplog << "Writing parameter file \"" << filename << "\"" << std::endl;
  prm.print_parameters(filename, ParameterHandler::PRM);
}

void
Sapphire::Utils::ParameterParser::write_template_parameters(
  const std::string &filename)
{
  LogStream::Prefix p("ParameterParser", saplog);
  saplog << "Writing template parameter file \"" << filename << "\""
         << std::endl;
  // The ParameterHandler is cleared, but the parameters stay the same, because
  // parse_aprameters() is not called.
  prm.clear();
  declare_parameters();
  prm.print_parameters(filename, ParameterHandler::PRM);

  // restore original state
  prm.parse_input(prm_file_name);
}

void
Sapphire::Utils::ParameterParser::declare_parameters()
{
  // VFP
  prm.enter_subsection("VFP");
  {
    prm.enter_subsection("Mesh");
    { // NOTE: This is a very strange syntax
      prm.declare_entry(
        "Point 1",
        "0., 0., 0.",
        dealii::Patterns::Anything(),
        "Two diagonally opposite corner points, Point 1 and  Point 2");
      prm.declare_entry(
        "Point 2",
        "1., 1., 1.",
        dealii::Patterns::Anything(),
        "Two diagonally opposite corner points Point 1 and  Point 2");
      prm.declare_entry("Number of cells",
                        "4, 4, 4",
                        dealii::Patterns::Anything(),
                        "Number of cells in each coordinate direction");
      prm.declare_entry(
        "Periodicity",
        "false, false, false",
        dealii::Patterns::Anything(),
        "Periodic boundaries in the three coordinate directions.");
    }
    prm.leave_subsection();

    prm.enter_subsection("Time stepping");
    {
      prm.declare_entry(
        "Method",
        "Crank-Nicolson",
        dealii::Patterns::Selection(
          "Forward Euler|Backward Euler|Crank-Nicolson|ERK4|LSERK"),
        "The time stepping method.");
      prm.declare_entry("Time step size",
                        "7.8125e-3",
                        dealii::Patterns::Double(),
                        "Duration of the simulation.");
      prm.declare_entry("Final time",
                        "0.4",
                        dealii::Patterns::Double(),
                        "Duration of the simulation.");
    }
    prm.leave_subsection();

    prm.enter_subsection("Expansion");
    {
      prm.declare_entry(
        "Expansion order",
        "0",
        dealii::Patterns::Integer(0),
        "The order of the expansion of the particle distribution function.");
    }
    prm.leave_subsection();

    prm.enter_subsection("Finite element");
    {
      prm.declare_entry("Polynomial degree",
                        "1",
                        dealii::Patterns::Integer(0),
                        "The degree of the shape functions (i.e. "
                        "the polynomials) of the finite element.");
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();

  saplog << "Declaring parameters" << std::endl;

  prm.enter_subsection("Output");
  {
    prm.declare_entry(
      "Results folder",
      "./results",
      dealii::Patterns::Anything(),
      "Path to the folder in which the simulation results will be stored. "
      "Without a trailing slash.");
    prm.declare_entry("Simulation identifier",
                      "",
                      dealii::Patterns::Anything(),
                      "Name of the simulation run. It will be "
                      "used to create a subdirectory "
                      "in the results folder.");
    prm.declare_entry("Base file name",
                      "solution",
                      dealii::Patterns::Anything(),
                      "The base file name for the output files.");
    prm.declare_entry("Number of digits for counter",
                      "4",
                      dealii::Patterns::Integer(0),
                      "The number of digits used for the counter in the "
                      "output file names.");
    prm.declare_entry(
      "Format",
      "vtu",
      dealii::Patterns::Selection("vtu|pvtu|hdf5"),
      "The format in which the simulation output will be stored.");
    prm.declare_entry("Output frequency",
                      "1",
                      dealii::Patterns::Integer(0),
                      "The frequence at which output files will "
                      "be written. (In units of time steps)");
  } // subsection Output
  prm.leave_subsection();

  prm.enter_subsection("Hydrodynamics");
  {
    prm.enter_subsection("Time stepping");
    {
      prm.declare_entry("Scheme",
                        "Forward Euler",
                        Patterns::Selection("Forward Euler|Explicit RK"),
                        "Time stepping scheme");
      prm.declare_entry("Time step",
                        "0.1",
                        Patterns::Double(0.0),
                        "Time step size");
      prm.declare_entry("End time",
                        "1.0",
                        Patterns::Double(0.0),
                        "End time of simulation");
    } // subsection Time stepping
    prm.leave_subsection();

    prm.enter_subsection("Finite element");
    {
      prm.declare_entry("Polynomial degree",
                        "1",
                        Patterns::Integer(1),
                        "Polynomial degree of finite element");
    } // subsection Finite element
    prm.leave_subsection();

    prm.enter_subsection("Mesh");
    {
      prm.declare_entry("Refinement level",
                        "1",
                        Patterns::Integer(0),
                        "Refinement level of mesh");
    } // subsection Mesh
    prm.leave_subsection();

    prm.enter_subsection("Numerical flux");
    {
      prm.declare_entry("Numerical flux",
                        "Lax-Friedrichs",
                        Patterns::Selection("Central|Upwind|Lax-Friedrichs"),
                        "Numerical flux");
    } // subsection Numerical flux
    prm.leave_subsection();

    prm.enter_subsection("Slope limiter");
    {
      prm.declare_entry("Slope limiter",
                        "No limiter",
                        Patterns::Selection(
                          "No limiter|Linear reconstruction|MinMod|MUSCL"),
                        "Slope limiter");
      prm.declare_entry(
        "Slope limiter criterion",
        "Never",
        Patterns::Selection("Never|Always|Generalized slope limiter"),
        "Criterion on which cells the slope limiter should apply");
    } // subsection Slope limiter
    prm.leave_subsection();

    prm.enter_subsection("Linear solver");
    {
      prm.declare_entry("Max iterations",
                        "1000",
                        Patterns::Integer(0),
                        "Maximum number of iterations");
      prm.declare_entry("Tolerance",
                        "1e-12",
                        Patterns::Double(0.0),
                        "Tolerance of solver");
    } // subsection Linear solver
    prm.leave_subsection();
  } // subsection Hydrodynamics
  prm.leave_subsection();
}

void
Sapphire::Utils::ParameterParser::parse_parameters()
{
  // VFP
  prm.enter_subsection("VFP");
  {
    prm.enter_subsection("Mesh");
    {
      // Two diagonally opposite corner points of the grid
      p1 = prm.get("Point 1");
      p2 = prm.get("Point 2");
      // Number of cells
      n_cells = prm.get("Number of cells");
      // Periodicity
      periodicity = prm.get("Periodicity");
    }
    prm.leave_subsection();

    prm.enter_subsection("Time stepping");
    {
      time_stepping_method = prm.get("Method");
      time_step            = prm.get_double("Time step size");
      final_time           = prm.get_double("Final time");
    }
    prm.leave_subsection();

    prm.enter_subsection("Expansion");
    {
      expansion_order = prm.get_integer("Expansion order");
    }
    prm.leave_subsection();

    prm.enter_subsection("Finite element");
    {
      polynomial_degree = prm.get_integer("Polynomial degree");
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();

  saplog << "Parsing parameters" << std::endl;

  std::string s;
  prm.enter_subsection("Output");
  {
    out_results_path         = prm.get("Results folder");
    out_simulation_id        = prm.get("Simulation identifier");
    out_base_file_name       = prm.get("Base file name");
    out_n_digits_for_counter = prm.get_integer("Number of digits for counter");

    std::string s;
    s = prm.get("Format");
    if (s == "vtu")
      out_format = Sapphire::Utils::OutputFormat::vtu;
    else if (s == "pvtu")
      out_format = Sapphire::Utils::OutputFormat::pvtu;
    else if (s == "hdf5")
      out_format = Sapphire::Utils::OutputFormat::hdf5;
    else
      AssertThrow(false, ExcNotImplemented());

    out_output_frequency = prm.get_integer("Output frequency");
  } // subsection Output
  prm.leave_subsection();

  prm.enter_subsection("Hydrodynamics");
  {
    prm.enter_subsection("Time stepping");
    {
      s = prm.get("Scheme");
      if (s == "Forward Euler")
        hdsolver_scheme = Sapphire::Hydro::TimeSteppingScheme::ForwardEuler;
      else if (s == "Explicit RK")
        hdsolver_scheme = Sapphire::Hydro::TimeSteppingScheme::ExplicitRK;
      else
        AssertThrow(false, ExcNotImplemented());

      hdsolver_time_step = prm.get_double("Time step");
      hdsolver_end_time  = prm.get_double("End time");
    } // subsection Time stepping
    prm.leave_subsection();

    prm.enter_subsection("Finite element");
    {
      hdsolver_fe_degree = prm.get_integer("Polynomial degree");
    } // subsection Finite element
    prm.leave_subsection();

    prm.enter_subsection("Mesh");
    {
      hdsolver_refinement_level = prm.get_integer("Refinement level");
    } // subsection Mesh
    prm.leave_subsection();

    prm.enter_subsection("Numerical flux");
    {
      s = prm.get("Numerical flux");
      if (s == "Central")
        hdsolver_flux_type = Sapphire::Hydro::FluxType::Central;
      else if (s == "Upwind")
        hdsolver_flux_type = Sapphire::Hydro::FluxType::Upwind;
      else if (s == "Lax-Friedrichs")
        hdsolver_flux_type = Sapphire::Hydro::FluxType::LaxFriedrichs;
      else
        AssertThrow(false, ExcNotImplemented());
    } // subsection Numerical flux
    prm.leave_subsection();

    prm.enter_subsection("Slope limiter");
    {
      s = prm.get("Slope limiter");
      if (s == "No limiter")
        hdsolver_limiter = Sapphire::Hydro::SlopeLimiter::NoLimiter;
      else if (s == "Linear reconstruction")
        hdsolver_limiter = Sapphire::Hydro::SlopeLimiter::LinearReconstruction;
      else if (s == "MinMod")
        hdsolver_limiter = Sapphire::Hydro::SlopeLimiter::MinMod;
      else if (s == "MUSCL")
        hdsolver_limiter = Sapphire::Hydro::SlopeLimiter::MUSCL;
      else
        AssertThrow(false, ExcNotImplemented());

      s = prm.get("Slope limiter criterion");
      if (s == "Never")
        hdsolver_limiter_criterion =
          Sapphire::Hydro::SlopeLimiterCriterion::Never;
      else if (s == "Always")
        hdsolver_limiter_criterion =
          Sapphire::Hydro::SlopeLimiterCriterion::Always;
      else if (s == "Generalized slope limiter")
        hdsolver_limiter_criterion =
          Sapphire::Hydro::SlopeLimiterCriterion::GerneralizedSlopeLimiter;
      else
        AssertThrow(false, ExcNotImplemented());

      // TODO: Do this here, or in init of HDSolverControl?
      if (hdsolver_limiter == Sapphire::Hydro::SlopeLimiter::NoLimiter)
        {
          Assert(hdsolver_limiter_criterion ==
                   Sapphire::Hydro::SlopeLimiterCriterion::Never,
                 ExcMessage("No slope limiter is selected, but a slope limiter "
                            "criterion is selected"));
        }
    } // subsection Slope limiter
    prm.leave_subsection();

    prm.enter_subsection("Linear solver");
    {
      hdsolver_max_iterations = prm.get_integer("Max iterations");
      hdsolver_tolerance      = prm.get_double("Tolerance");
    } // subsection Linear solver
    prm.leave_subsection();
  } // subsection Hydrodynamics
  prm.leave_subsection();
}
