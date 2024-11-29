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
     * Like
     * @dealref{InterpolatedUniformGridData,classFunctions_1_1InterpolatedUniformGridData}
     * but for with the option to update data.
     *
     * @todo Implement in @dealii.
     * @tparam dim Dimension of the data in the file.
     */
    template <int dim>
    class InterpolatedUniformGridData2 : public Function<dim>
    {
    public:
      /**
       * Constructor
       * @param interval_endpoints The left and right end points of the
       * (uniformly subdivided) intervals in each of the coordinate directions.
       * @param n_subintervals The number of subintervals in each coordinate
       * direction. A value of one for a coordinate means that the interval is
       * considered as one subinterval consisting of the entire range. A value
       * of two means that there are two subintervals each with one half of the
       * range, etc.
       * @param data_values A dim-dimensional table of data at each of the mesh
       * points defined by the coordinate arrays above. Note that the Table
       * class has a number of conversion constructors that allow converting
       * other data types into a table where you specify this argument.
       */
      InterpolatedUniformGridData2(
        const std::array<std::pair<double, double>, dim> &interval_endpoints,
        const std::array<unsigned int, dim>              &n_subintervals,
        const dealii::Table<dim, double>                 &data_values);

      /**
       * Like the previous constructor, but take the arguments as rvalue
       * references and *move*, instead of *copy* the data. This is often useful
       * in cases where the data stored in these tables is large and the
       * information used to initialize the current object is no longer needed
       * separately. In other words, there is no need to keep the original
       * object from which this object could copy its information, but it might
       * as well take over ("move") the data.
       *
       * Moving data also enables using tables that are located in shared memory
       * between multiple MPI processes, rather than copying the data from
       * shared memory into local memory whenever one creates an
       * InterpolatedUniformGridData object. See the
       * TableBase::replicate_across_communicator() function on how to share a
       * data set between multiple processes.
       */
      InterpolatedUniformGridData2(
        std::array<std::pair<double, double>, dim> &&interval_endpoints,
        std::array<unsigned int, dim>              &&n_subintervals,
        dealii::Table<dim, double>                 &&data_values);

      /**
       * Empty constructor
       */
      InterpolatedUniformGridData2() = default;

      /**
       * Compute the value of the function set by bilinear interpolation of the
       * given data set.
       *
       * @param p The point at which the function is to be evaluated.
       * @param component The vector component. Since this function is scalar,
       * only zero is a valid argument here.
       * @return The interpolated value at this point. If the point lies outside
       * the set of coordinates, the function is extended by a constant.
       */
      virtual double
      value(const Point<dim>  &p,
            const unsigned int component = 0) const override;

      /**
       * Compute the gradient of the function set by bilinear interpolation of
       * the given data set.
       *
       * @param p The point at which the function is to be evaluated.
       * @param component The vector component. Since this function is scalar,
       *   only zero is a valid argument here.
       * @return The gradient of the interpolated function at this point. If the
       *   point lies outside the set of coordinates, the function is extended
       *   by a constant whose gradient is then of course zero.
       */
      virtual Tensor<1, dim>
      gradient(const Point<dim>  &p,
               const unsigned int component = 0) const override;

      /**
       * Return an estimate for the memory consumption, in bytes, of this
       * object.
       */
      virtual std::size_t
      memory_consumption() const override;

      /**
       * Return a reference to the internally stored data.
       */
      const Table<dim, double> &
      get_data() const;

      /**
       * Updates the data used by the function.
       * @param interval_endpoints The left and right end points of the
       * (uniformly subdivided) intervals in each of the coordinate directions.
       * @param n_subintervals The number of subintervals in each coordinate
       * direction. A value of one for a coordinate means that the interval is
       * considered as one subinterval consisting of the entire range. A value
       * of two means that there are two subintervals each with one half of the
       * range, etc.
       * @param data_values A dim-dimensional table of data at each of the mesh
       * points defined by the coordinate arrays above. Note that the Table
       * class has a number of conversion constructors that allow converting
       * other data types into a table where you specify this argument.
       */
      void
      set_data(
        const std::array<std::pair<double, double>, dim> &interval_endpoints,
        const std::array<unsigned int, dim>              &n_subintervals,
        const dealii::Table<dim, double>                 &data_values);

      /**
       * Like the previous function, but take the arguments as rvalue
       * references and *move*, instead of *copy* the data. This is often useful
       * in cases where the data stored in these tables is large and the
       * information used to initialize the current object is no longer needed
       * separately. In other words, there is no need to keep the original
       * object from which this object could copy its information, but it might
       * as well take over ("move") the data.
       *
       * Moving data also enables using tables that are located in shared memory
       * between multiple MPI processes, rather than copying the data from
       * shared memory into local memory whenever one creates an
       * InterpolatedUniformGridData object. See the
       * TableBase::replicate_across_communicator() function on how to share a
       * data set between multiple processes.
       */
      void
      set_data(std::array<std::pair<double, double>, dim> &&interval_endpoints,
               std::array<unsigned int, dim>              &&n_subintervals,
               dealii::Table<dim, double>                 &&data_values);

    private:
      /**
       * The set of interval endpoints in each of the coordinate directions.
       */
      std::array<std::pair<double, double>, dim> interval_endpoints;

      /**
       * The number of subintervals in each of the coordinate directions.
       */
      std::array<unsigned int, dim> n_subintervals;

      /**
       * The data that is to be interpolated.
       */
      Table<dim, double> data_values;
    };



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
