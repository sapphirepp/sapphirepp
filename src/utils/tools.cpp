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
 * @file tools.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement functions in namespace @ref sapphirepp::Utils::Tools.
 */


#include "tools.h"



std::vector<double>
sapphirepp::Utils::Tools::create_linear_range(const double       start,
                                              const double       stop,
                                              const unsigned int num)
{
  std::vector<double> values(num);
  if (num == 0)
    return values;

  for (unsigned int i = 0; i < num; ++i)
    values[i] = start + i * (stop - start) / (num - 1);
  // Sanitize
  values[0]       = start;
  values[num - 1] = stop;
  return values;
}
