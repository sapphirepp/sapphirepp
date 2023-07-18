#include "parameter-parser.h"

#include <iostream>

#include "sapphire-logstream.h"

Sapphire::Utils::ParameterParser::ParameterParser(
  const std::string &prm_file_name)
{
  LogStream::Prefix p("ParameterParser", saplog);
  declare_parameters();
  saplog << "Parsing input file \"" << prm_file_name << "\"" << std::endl;
  prm.parse_input(prm_file_name);
  parse_parameters();
  prm.log_parameters(saplog);
}

void
Sapphire::Utils::ParameterParser::declare_parameters()
{
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
}

void
Sapphire::Utils::ParameterParser::parse_parameters()
{
  saplog << "Parsing parameters" << std::endl;
  std::string s;

  prm.enter_subsection("Output");
  {
    out_results_path         = prm.get("Results folder");
    out_simulation_id        = prm.get("Simulation identifier");
    out_base_file_name       = prm.get("Base file name");
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
}