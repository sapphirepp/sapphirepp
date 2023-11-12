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
 * @file sapphire-logstream.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define saplog
 * @version 0.1
 * @date 2023-07-17
 */

#ifndef UTILS_SAPPHIREPPLOGSTREAM_H
#define UTILS_SAPPHIREPPLOGSTREAM_H

#include <deal.II/base/logstream.h>

#include <string>

namespace Sapphire
{
  namespace Utils
  {

    class SapphireppLogStream : public dealii::LogStream
    {
    public:
      SapphireppLogStream();
    };

  } // namespace Utils

  extern Sapphire::Utils::SapphireppLogStream saplog;
} // namespace Sapphire
#endif
