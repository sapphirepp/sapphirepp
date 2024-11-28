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
 * @file grid-data-function.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::Utils::GridDataFunction
 */

#ifndef UTILS_GRIDDATAFUNCTION_H
#define UTILS_GRIDDATAFUNCTION_H

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/tensor.h>

#include <filesystem>
#include <string>
#include <vector>

namespace sapphirepp
{
  namespace Utils
  {
    using namespace dealii;

    /**
     * @brief @dealref{Function} created from data in a file.
     *
     * This function reads data from a file and uses the
     * @dealref{InterpolatedTensorProductGridData,classFunctions_1_1InterpolatedTensorProductGridData}
     * class to create a function with multiple components.
     *
     * @tparam dim Dimension of the data in the file.
     * @tparam spacedim Space dimension of the function.
     */

    /**
     * @brief Constructor
     */
    template <unsigned int dim, unsigned int spacedim = dim>
    class GridDataFunction : public Function<spacedim>
    {
    public:
      /** @brief Constructor */
      GridDataFunction(const unsigned int n_components = 1,
                       const double       inital_time  = 0.0);



      virtual double
      value(const Point<spacedim> &p,
            const unsigned int     component = 0) const override;



      virtual Tensor<1, spacedim>
      gradient(const Point<spacedim> &p,
               const unsigned int     component = 0) const override;



      virtual void
      set_time(const double new_time) override;
    };

  } // namespace Utils
} // namespace sapphirepp
#endif
