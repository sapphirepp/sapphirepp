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
 * @file tools.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define functions in namespace @ref sapphirepp::Utils::Tools.
 */

#ifndef UTILS_TOOLS_H
#define UTILS_TOOLS_H


#include <vector>



namespace sapphirepp
{
  namespace Utils
  {

    /**
     * @namespace sapphirepp::Utils::Tools
     * @brief Namespace for general utility functions
     */
    namespace Tools
    {



      /**
       * @brief Create a vector with `num` linear spaced values.
       *
       * @param start Starting value
       * @param stop End value
       * @param num Number of values
       * @return std::vector<double>
       */
      std::vector<double>
      create_linear_range(const double       start,
                          const double       stop,
                          const unsigned int num);

    } // namespace Tools
  }   // namespace Utils
} // namespace sapphirepp
#endif
