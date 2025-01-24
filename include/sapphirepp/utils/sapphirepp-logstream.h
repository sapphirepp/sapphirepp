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
 * @file sapphirepp-logstream.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define  @ref sapphirepp::Utils::SapphireppLogStream and
 *        @ref sapphirepp::saplog
 */

#ifndef UTILS_SAPPHIREPPLOGSTREAM_H
#define UTILS_SAPPHIREPPLOGSTREAM_H

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>

#include <exception>
#include <fstream>
#include <string>

/**
 * @namespace sapphirepp
 * @brief Namespace for @sapphire
 */
namespace sapphirepp
{

  /**
   * @namespace sapphirepp::Utils
   * @brief Namespace for utility functions
   *
   * This namespace contains utility functions for @sapphire, that are not
   * directly related to the VFP module.
   */
  namespace Utils
  {

    /**
     * @brief @dealref{LogStream} for @sapphire
     *
     * The Stream is prefixed with `sapphirepp`.
     * In case of parallel execution, all but the first process are prefixed
     * with `mpi<rank>`.
     *
     * For more details see @dealref{LogStream} documentation
     */
    class SapphireppLogStream : public dealii::LogStream
    {
    public:
      /** @brief Constructor */
      SapphireppLogStream();


      /** @brief Destructor */
      ~SapphireppLogStream();


      /**
       * @brief Initialize the log stream
       *
       * This function must only be called after MPI initialization.
       * It prefixes the log stream with `MPI<rank>` for all but the first
       * rank. In addition, it shows a start-up message.
       *
       * @param depth_console The verbosity of the console output
       *        - `0` silence the program
       *        - `1` shows only the start-up message
       *        - `2` show progress
       *        - `>2` show different levels of debug messages
       * @param enable_mpi_output If `true`, all mpi processes will output to
       *        the console, otherwise only the first process will show messages
       */
      void
      init(const unsigned int depth_console     = 2,
           const bool         enable_mpi_output = false);


      /**
       * @brief Initialize the log stream with command line options
       *
       * This function must only be called after MPI initialization.
       * To show all options use the `--help/-h` flag.
       *
       * @param argc Number of commandline arguments
       * @param argv Commandline arguments
       */
      void
      init(const int argc, const char *const argv[]);


      /**
       * @brief Attach a logfile to the stream
       *
       * @param filepath Path to the logfile
       * @param depth_file The verbosity of the file output
       * @param enable_mpi_output Enable output of all mpi processes
       */
      void
      attach_file(const std::string &filepath,
                  const unsigned int depth_file,
                  const bool         enable_mpi_output = false);


      /**
       * @brief Get the verbosity of the console output
       *
       * @return unsigned int verbosity
       */
      unsigned int
      get_verbosity();


      /**
       * @brief Prints an error message to `std:err` and this LogStream.
       *
       * This has the advantage, that error messages are also saved in the
       * logfile.
       *
       * @param exc Exception
       */
      void
      print_error(const std::exception &exc);



    private:
      /** Output stream to logfile */
      std::ofstream log_file;
    };

  } // namespace Utils

  /**
   * @brief The standard @ref Utils::SapphireppLogStream "log stream"
   *        for @sapphire
   */
  extern sapphirepp::Utils::SapphireppLogStream saplog;
} // namespace sapphirepp
#endif
