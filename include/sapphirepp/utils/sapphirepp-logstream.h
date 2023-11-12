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

/// @file sapphirepp-logstream.h
/// @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
/// @brief Define saplog

#ifndef UTILS_SAPPHIREPPLOGSTREAM_H
#define UTILS_SAPPHIREPPLOGSTREAM_H

#include <deal.II/base/logstream.h>

namespace Sapphire
{
  namespace Utils
  {

    /// @brief LogStream for @sapphire
    ///
    /// The Stream is prefixed with `Sapphire`.
    /// In case of parallel execution, all but the first process are prefixed
    /// with `mpi<rank>`.
    ///
    /// For more details see @dealref{LogStream} documentation
    class SapphireppLogStream : public dealii::LogStream
    {
    public:
      /// @brief Constructor
      SapphireppLogStream();

      /// @brief Initialize the log stream
      /// This function must only be called after MPI initialization.
      /// It prefixes the log stream with `MPI<rank>` for all but the first
      /// rank. In addition, it shows a start-up message.
      void
      init();
    };

  } // namespace Utils

  /// @brief The standard log stream for @sapphire
  extern Sapphire::Utils::SapphireppLogStream saplog;
} // namespace Sapphire
#endif
