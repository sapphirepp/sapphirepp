#include "vfp-solver-control.h"

VFPEquation::VFPSolverControl::VFPSolverControl(const std::string& file_path)
    : parameter_file{file_path} {
  declare_parameters();
  parse_parameters();
  get_parameters();
}

void VFPEquation::VFPSolverControl::print_settings(std::ostream& os) const {
  os << "Compile time parameters: " << std::endl;
  os << "	Dimension Configuration Space: "
     << VFPSolverControl::dim_configuration_space << "\n";
  os << "	" << terms << std::endl;
  os << "Runtime parameters:" << std::endl;
  parameter_handler.print_parameters(os, dealii::ParameterHandler::ShortText);
  os << std::endl;
}

void VFPEquation::VFPSolverControl::declare_parameters() {
  parameter_handler.enter_subsection("Mesh");
  {  // NOTE: This is a very strange syntax
    parameter_handler.declare_entry("Number of refinements", "6",
                                    dealii::Patterns::Integer(0),
                                    "Number of global mesh refinement steps");
  }
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection("Time stepping");
  {
    parameter_handler.declare_entry("Time step size", "7.8125e-3",
                                    dealii::Patterns::Double(),
                                    "Duration of the simulation.");
    parameter_handler.declare_entry("Final time", "0.4",
                                    dealii::Patterns::Double(),
                                    "Duration of the simulation.");
  }
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection("Expansion");
  {
    parameter_handler.declare_entry(
        "Expansion order", "0", dealii::Patterns::Integer(0),
        "The order of the expansion of the particle distribution function.");
  }
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection("Finite element");
  {
    parameter_handler.declare_entry("Polynomial degree", "1",
                                    dealii::Patterns::Integer(0),
                                    "The degree of the shape functions (i.e. "
                                    "the polynomials) of the finite element.");
  }
  parameter_handler.leave_subsection();
}

void VFPEquation::VFPSolverControl::parse_parameters() {
  parameter_handler.parse_input(parameter_file);
}

void VFPEquation::VFPSolverControl::get_parameters() {
  parameter_handler.enter_subsection("Mesh");
  { num_refinements = parameter_handler.get_integer("Number of refinements"); }
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection("Time stepping");
  {
    time_step = parameter_handler.get_double("Time step size");
    final_time = parameter_handler.get_double("Final time");
  }
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection("Expansion");
  { expansion_order = parameter_handler.get_integer("Expansion order"); }
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection("Finite element");
  { polynomial_degree = parameter_handler.get_integer("Polynomial degree"); }
  parameter_handler.leave_subsection();
}