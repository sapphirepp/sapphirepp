#include "parameter_parser.h"

#include <iostream>

Sapphire::Utils::ParameterParser::ParameterParser(
    const std::string &prm_file_name) {
  declare_parameters();
  prm.parse_input(prm_file_name);
  prm.print_parameters(std::cout, ParameterHandler::PRM);
  parse_parameters();
}

void Sapphire::Utils::ParameterParser::declare_parameters() {
  prm.enter_subsection("Output");
  {
    prm.declare_entry(
        "Results folder", "./results", dealii::Patterns::Anything(),
        "Path to the folder in which the simulation results will be stored. "
        "Without a trailing slash.");
    prm.declare_entry("Simulation identifier", "", dealii::Patterns::Anything(),
                      "Name of the simulation run. It will be "
                      "used to create a subdirectory "
                      "in the results folder.");
    prm.declare_entry("Base file name", "solution",
                      dealii::Patterns::Anything(),
                      "The base file name for the output files.");
    prm.declare_entry("Number of digits for counter", "4",
                      dealii::Patterns::Integer(0),
                      "The number of digits used for the counter in the "
                      "output file names.");
    prm.declare_entry(
        "Format", "vtu", dealii::Patterns::Selection("vtu|pvtu|hdf5"),
        "The format in which the simulation output will be stored.");
    prm.declare_entry("Output frequency", "1", dealii::Patterns::Integer(0),
                      "The frequence at which output files will "
                      "be written. (In units of time steps)");
  } // subsection Output
  prm.leave_subsection();

  prm.enter_subsection("Hydrodynamics");
  {

    prm.enter_subsection("Time stepping");
    {
      prm.declare_entry("Scheme", "Forward Euler",
                        Patterns::Selection("Forward Euler|Explicit RK"),
                        "Time stepping scheme");
      prm.declare_entry("Time step", "0.1", Patterns::Double(0.0),
                        "Time step size");
      prm.declare_entry("End time", "1.0", Patterns::Double(0.0),
                        "End time of simulation");
    } // subsection Time stepping
    prm.leave_subsection();

    prm.enter_subsection("Finite element");
    {
      prm.declare_entry("Polynomial degree", "1", Patterns::Integer(1),
                        "Polynomial degree of finite element");
    } // subsection Finite element
    prm.leave_subsection();

    prm.enter_subsection("Mesh");
    {
      prm.declare_entry("Refinement level", "1", Patterns::Integer(0),
                        "Refinement level of mesh");
    } // subsection Mesh
    prm.leave_subsection();

    prm.enter_subsection("Numerical flux");
    {
      prm.declare_entry("Numerical flux", "Lax-Friedrichs",
                        Patterns::Selection("Central|Upwind|Lax-Friedrichs"),
                        "Numerical flux");
    } // subsection Numerical flux
    prm.leave_subsection();

    prm.enter_subsection("Slope limiter");
    {
      prm.declare_entry(
          "Slope limiter", "No limiter",
          Patterns::Selection("No limiter|Linear reconstruction|MinMod|MUSCL"),
          "Slope limiter");
      prm.declare_entry(
          "Slope limiter criterion", "Never",
          Patterns::Selection("Never|Always|Generalized slope limiter"),
          "Criterion on which cells the slope limiter should apply");
    } // subsection Slope limiter
    prm.leave_subsection();

    prm.enter_subsection("Linear solver");
    {
      prm.declare_entry("Max iterations", "1000", Patterns::Integer(0),
                        "Maximum number of iterations");
      prm.declare_entry("Tolerance", "1e-12", Patterns::Double(0.0),
                        "Tolerance of solver");
    } // subsection Linear solver
    prm.leave_subsection();
  } // subsection Hydrodynamics
  prm.leave_subsection();
}

void Sapphire::Utils::ParameterParser::parse_parameters() {
  std::string s;

  prm.enter_subsection("Output");
  {
    out_results_path = prm.get("Results folder");
    out_simulation_id = prm.get("Simulation identifier");
    out_base_file_name = prm.get("Base file name");
    out_n_digits_for_counter = prm.get_integer("Number of digits for counter");

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
      hdsolver_end_time = prm.get_double("End time");
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
      if (hdsolver_limiter == Sapphire::Hydro::SlopeLimiter::NoLimiter) {
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
      hdsolver_tolerance = prm.get_double("Tolerance");
    } // subsection Linear solver
    prm.leave_subsection();
  } // subsection Hydrodynamics
  prm.leave_subsection();
}