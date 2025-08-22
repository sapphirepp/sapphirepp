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
 * @file sapphirepp-logstream.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement @ref sapphirepp::Utils::SapphireppLogStream
 */

#include "sapphirepp-logstream.h"

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <boost/program_options.hpp>

#include <petscsys.h>

#include <mpi.h>

#include <iostream>
#include <limits>

#include "version.h"



sapphirepp::Utils::SapphireppLogStream sapphirepp::saplog;



sapphirepp::Utils::SapphireppLogStream::SapphireppLogStream()
  : dealii::LogStream()
{
  this->pop();
  this->push("Sapphire");
}



sapphirepp::Utils::SapphireppLogStream::~SapphireppLogStream()
{
  // Ensure that the log_file is closed.
  // I think we do not need the same precautions as for the deallog:
  // https://www.dealii.org/current/doxygen/deal.II/logstream_8cc_source.html#l00088
  if (log_file.is_open())
    {
      if (dealii::deallog.has_file())
        dealii::deallog.detach();
      this->detach();
      log_file.close();
    }
}



void
sapphirepp::Utils::SapphireppLogStream::init(const unsigned int depth_console,
                                             const bool enable_mpi_output)
{
  int is_initialized;
  MPI_Initialized(&is_initialized);
  Assert(is_initialized,
         dealii::ExcMessage("saplog.init() must be called after MPI "
                            "initialization!"));
  this->depth_console(depth_console);
  this->pop();
  const unsigned int mpi_rank =
    dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (mpi_rank > 0)
    {
      this->push("MPI" + dealii::Utilities::to_string(mpi_rank, 3));
      if (!enable_mpi_output)
        this->depth_console(0);
    }
  this->push("Sapphire");

  *this << "Start Sapphire++ v" << SAPPHIREPP_VERSION << " with "
        << dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
        << " MPI process(es) [" << dealii::Utilities::System::get_date() << " "
        << dealii::Utilities::System::get_time() << "]" << std::endl;
}



void
sapphirepp::Utils::SapphireppLogStream::init(const int         argc,
                                             const char *const argv[])
{
  int is_initialized;
  MPI_Initialized(&is_initialized);
  Assert(is_initialized,
         dealii::ExcMessage("saplog.init() must be called after MPI "
                            "initialization!"));

  namespace po = boost::program_options;

  po::options_description desc("saplog options");
  desc.add_options()("help,h", "produce help message")(
    "verbosity,v",
    po::value<unsigned int>()->default_value(2),
    "verbosity of the console output\n"
    "  0 - silence the program\n"
    "  1 - show only start-up message\n"
    "  2 - show progress\n"
    "  >2 - show different levels of debug messages")(
    "verbose-mpi", "show output from all mpi processes")(
    "logfile,l", po::value<std::string>(), "output to logfile")(
    "logfile-verbosity",
    po::value<unsigned int>()->default_value(
      std::numeric_limits<unsigned int>::max()),
    "verbosity of the logfile output")("deallog",
                                       po::value<unsigned int>()->default_value(
                                         0),
                                       "verbosity of deallog");

  po::variables_map vm;
  auto              parsed = po::command_line_parser(argc, argv)
                  .options(desc)
                  .allow_unregistered()
                  .run();
  po::store(parsed, vm);
  po::notify(vm);

  if (vm.count("help"))
    {
      std::cout << std::endl
                << desc << std::endl
                << "  -help \t shows all MPI/PETSc options" << std::endl
                << std::endl;
    }

  // RemoveOptions from PETSC options to avoid warnings
  PetscOptionsClearValue(nullptr, "-h");
  PetscOptionsClearValue(nullptr, "--help");
  PetscOptionsClearValue(nullptr, "-v");
  PetscOptionsClearValue(nullptr, "--verbosity");
  PetscOptionsClearValue(nullptr, "--verbose-mpi");
  PetscOptionsClearValue(nullptr, "-l");
  PetscOptionsClearValue(nullptr, "--logfile");
  PetscOptionsClearValue(nullptr, "--logfile-verbosity");
  PetscOptionsClearValue(nullptr, "--deallog");

  if (vm.count("logfile"))
    {
      this->attach_file(vm["logfile"].as<std::string>(),
                        vm["logfile-verbosity"].as<unsigned int>(),
                        vm.count("verbose-mpi"));
    }
  this->init(vm["verbosity"].as<unsigned int>(), vm.count("verbose-mpi"));

  dealii::deallog.depth_console(vm["deallog"].as<unsigned int>());
  if (this->has_file())
    {
      dealii::deallog.attach(log_file, false);
      dealii::deallog.depth_file(vm["deallog"].as<unsigned int>());
    }
}



void
sapphirepp::Utils::SapphireppLogStream::attach_file(
  const std::string &filepath,
  const unsigned int depth_file,
  const bool         enable_mpi_output)
{
  Assert(!log_file.is_open(),
         dealii::ExcMessage("You try to attach a log_file, "
                            "but the log_filestream is already open. "
                            "You can attach only one log_file!"));
  int is_initialized;
  MPI_Initialized(&is_initialized);
  Assert(is_initialized,
         dealii::ExcMessage("saplog.init() must be called after MPI "
                            "initialization!"));

  const unsigned int mpi_rank =
    dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if ((mpi_rank > 0) && !enable_mpi_output)
    return;

  log_file.open(filepath, std::ios::app);
  AssertThrow(log_file.is_open(), dealii::ExcFileNotOpen(filepath));
  saplog.attach(log_file, mpi_rank == 0);
  this->depth_file(depth_file);
}



void
sapphirepp::Utils::SapphireppLogStream::print_error(const std::exception &exc)
{
  std::cerr << std::endl;
  std::cerr << "\n"
            << "----------------------------------------------------" << "\n"
            << "Exception on processing: " << "\n"
            << exc.what() << "\n"
            << "Aborting!" << "\n"
            << "----------------------------------------------------"
            << std::endl;

  *this << std::endl;
  // Always print the error to the log file
  const unsigned int depth_file =
    this->depth_file(std::numeric_limits<unsigned int>::max());
  *this << "\n"
        << "----------------------------------------------------" << "\n"
        << "Exception on processing: " << "\n"
        << exc.what() << "\n"
        << "Aborting!" << "\n"
        << "----------------------------------------------------" << std::endl;
  this->depth_file(depth_file);
}



void
sapphirepp::Utils::SapphireppLogStream::print_warning(
  const std::string &warning)
{
  if (this->get_verbosity_console() > 0)
    {
      std::cerr << std::endl;
      std::cerr << "WARNING: " << warning << std::endl;
    }

  *this << std::endl;
  // Always print the warning to the log file
  const unsigned int depth_file =
    this->depth_file(std::numeric_limits<unsigned int>::max());
  *this << "\n"
        << "WARNING: " << warning << std::endl;
  this->depth_file(depth_file);
}



unsigned int
sapphirepp::Utils::SapphireppLogStream::get_verbosity_console()
{
  unsigned int n = this->depth_console(0);
  this->depth_console(n);
  return n;
}



unsigned int
sapphirepp::Utils::SapphireppLogStream::get_verbosity_file()
{
  unsigned int n = this->depth_file(0);
  this->depth_file(n);
  return n;
}



std::ostream &
sapphirepp::Utils::SapphireppLogStream::to_ostream(const unsigned int verbosity)
{
  if (this->has_file() && (this->get_verbosity_file() >= verbosity))
    return this->get_file_stream();

  if (this->get_verbosity_console() >= verbosity)
    return this->get_console();

  // If the logfile is activated, we use write to the logfile,
  // otherwise the stream is discarded.
  return log_file;
}



dealii::ConditionalOStream
sapphirepp::Utils::SapphireppLogStream::to_condition_ostream(
  const unsigned int verbosity)
{
  if (this->has_file())
    return dealii::ConditionalOStream(this->get_file_stream(),
                                      this->get_verbosity_file() >= verbosity);

  return dealii::ConditionalOStream(this->get_console(),
                                    this->get_verbosity_console() >= verbosity);
}



sapphirepp::Utils::SapphireppLogStream::operator std::ostream &()
{
  return this->to_ostream();
}



sapphirepp::Utils::SapphireppLogStream::operator dealii::ConditionalOStream()
{
  return this->to_condition_ostream();
}
