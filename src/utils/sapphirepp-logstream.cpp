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

#include <string>

#include "version.h"



sapphirepp::Utils::SapphireppLogStream sapphirepp::saplog;



sapphirepp::Utils::SapphireppLogStream::SapphireppLogStream()
  : dealii::LogStream()
{
  this->pop();
  this->push("Sapphire");
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
      if (enable_mpi_output == false)
        this->depth_console(0);
    }
  this->push("Sapphire");

  *this << "Start Sapphire++ v" << SAPPHIREPP_VERSION << " with "
        << dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
        << " MPI process(es) [" << dealii::Utilities::System::get_date() << " "
        << dealii::Utilities::System::get_time() << "]" << std::endl;
}



void
sapphirepp::Utils::SapphireppLogStream::init(int argc, char *argv[])
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
    "verbose-mpi", "show output from all mpi processes");

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

  this->init(vm["verbosity"].as<unsigned int>(), vm.count("verbose-mpi"));
}



unsigned int
sapphirepp::Utils::SapphireppLogStream::get_verbosity()
{
  unsigned int n = this->depth_console(0);
  this->depth_console(n);
  return n;
}
