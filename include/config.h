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
 * @file config.h
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define and implement user defined physical setup
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <cmath>
#include <vector>

#include "pde-system.h"
#include "sapphirepp-logstream.h"
#include "vfp-flags.h"



namespace sapphirepp
{
  /**
   * @brief Class to implement user defined parameters
   */
  class PhysicalParameters
  {
  public:
    /** [Define runtime parameter] */
    // !!!EDIT HERE!!!
    /** [Define runtime parameter] */



    /** Constructor */
    PhysicalParameters() = default;



    /**
     * @brief Declare parameters in parameter file
     *
     * @param prm @dealref{ParameterHandler}
     */
    void
    declare_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix pre1("Startup", saplog);
      dealii::LogStream::Prefix pre2("Physical parameters", saplog);
      saplog << "Declaring parameters" << std::endl;
      prm.enter_subsection("Physical parameters");

      /** [Declare runtime parameter] */
      // !!!EDIT HERE!!!
      /** [Declare runtime parameter] */

      prm.leave_subsection();
    }



    /**
     * @brief Parse parameters from parameter file
     *
     * @param prm @dealref{ParameterHandler}
     */
    void
    parse_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix pre1("Startup", saplog);
      dealii::LogStream::Prefix pre2("PhysicalParameters", saplog);
      saplog << "Parsing parameters" << std::endl;
      prm.enter_subsection("Physical parameters");

      /** [Parse runtime parameter] */
      // !!!EDIT HERE!!!
      /** [Parse runtime parameter] */

      prm.leave_subsection();
    }
  };



  namespace VFP
  {
    /** [Dimension] */
    // !!!EDIT HERE!!!
    /** Specify reduced phase space dimension \f$ (\mathbf{x}, p) \f$ */
    static constexpr unsigned int dimension = 2;
    /** [Dimension] */



    /** [VFP Flags] */
    // !!!EDIT HERE!!!
    /** Specify which terms of the VFP equation should be active */
    static constexpr VFPFlags vfp_flags = VFPFlags::spatial_advection |
                                          VFPFlags::time_independent_fields |
                                          VFPFlags::source;
    /** [VFP Flags] */



    /**
     * @brief Initial condition
     *
     * @tparam dim Dimension of the reduced phase space \f$ (\mathbf{x}, p) \f$
     */
    template <unsigned int dim>
    class InitialValueFunction : public dealii::Function<dim>
    {
    public:
      /**
       * @brief Constructor
       *
       * @param physical_parameters User defined runtime parameters
       * @param exp_order Maximum expansion order \f$ l_{\rm max} \f$
       */
      InitialValueFunction(const PhysicalParameters &physical_parameters,
                           const unsigned int        exp_order)
        : dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
        , prm{physical_parameters}
        , lms_indices{PDESystem::create_lms_indices(exp_order)}
      {}



      /**
       * @brief Evaluate the initial condition at point `p`
       *
       * Return a vector of values corresponding to the initial condition of
       * the coefficients \f$ f_{i(l,m,s)}(t=0, \mathbf{x}, p) \f$ of the
       * expansion of the distribution function in spherical harmonics.
       *
       * The zeroth component of values corresponds to \f$ f_{000} \f$, the
       * first component to \f$ f_{110} \f$ etc. The mapping between the
       * system index \f$ i \f$ and the spherical harmonic indices \f$
       * (l,m,s) \f$ is given by `lms_indices`:
       *
       * ```cpp
       *  const unsigned int l = lms_indices[i][0];
       *  const unsigned int m = lms_indices[i][1];
       *  const unsigned int s = lms_indices[i][2];
       * ```
       *
       * The point in reduced phase space is given by `p`. The first
       * `dim_cs` components correspond to the spatial coordinates \f$
       * \mathbf{x} \f$ and if momentum term is activated the last component
       * to the momentum coordinate \f$ p \f$:
       *
       * ```cpp
       *  const double x = p[0];
       *  const double y = p[1];
       *  ...
       *  const double p = p[dim-1];
       * ```
       *
       * @param p Point in reduced phase space
       * @param f Return vector \f$ f_{i(l,m,s)}(t=0, \mathbf{x}, p) \f$
       * @see @dealref{Function::vector_value(),classFunction,ae316ebc05d21989d573024f8a23c49cb}
       */
      void
      vector_value(const dealii::Point<dim> &p,
                   dealii::Vector<double>   &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(p); // suppress compiler warning

        for (unsigned int i = 0; i < f.size(); ++i)
          {
            /** [Initial value] */
            // !!!EDIT HERE!!!
            f[i] = 0.;
            /** [Initial value] */
          }
      }



    private:
      /** User defined runtime parameters */
      const PhysicalParameters prm;
      /**
       * Map between system index \f$ i \f$ and spherical harmonic indices
       * \f$ (l,m,n) \f$
       */
      const std::vector<std::array<unsigned int, 3>> lms_indices;
    };



    /**
     * @brief Scattering frequency
     *
     * @tparam dim Dimension of the reduced phase space \f$ (\mathbf{x}, p) \f$
     */
    template <unsigned int dim>
    class ScatteringFrequency : public dealii::Function<dim>
    {
    public:
      /**
       * @brief Constructor
       *
       * @param physical_parameters User defined runtime parameters
       */
      ScatteringFrequency(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(1)
        , prm{physical_parameters}
      {}



      /**
       * @brief Evaluate the scattering frequency
       *
       * Return a vector of values corresponding to the scattering frequency
       * \f$ \nu(\mathbf{x}, p) \f$ at all `points` in reduced phase
       * space. `points` is a vector of points in reduced phase space:
       *
       * ```cpp
       * for (unsigned int i = 0; i < points.size(); ++i)
       *  {
       *    const double x = points[i][0];
       *    const double y = points[i][1];
       *    ...
       *    const double p = points[i][dim-1];
       *  }
       * ```
       *
       * @param points List of points in reduced phase space
       * @param scattering_frequencies List of scattering frequencies
       * @param component Unused parameter
       * @see @dealref{Function::value_list(),classFunction,abe86ee7f7f12cf4041d1e714c0fb42f3}
       */
      void
      value_list(const std::vector<dealii::Point<dim>> &points,
                 std::vector<double>                   &scattering_frequencies,
                 const unsigned int component = 0) const override
      {
        AssertDimension(scattering_frequencies.size(), points.size());
        static_cast<void>(component); // suppress compiler warning

        for (unsigned int i = 0; i < points.size(); ++i)
          {
            /** [Scattering frequency] */
            // !!!EDIT HERE!!!
            scattering_frequencies[i] = 0.;
            /** [Scattering frequency] */
          }
      }



    private:
      /** User defined runtime parameters */
      const PhysicalParameters prm;
    };



    /**
     * @brief Source term
     *
     * @tparam dim Dimension of the reduced phase space \f$ (\mathbf{x}, p) \f$
     */
    template <unsigned int dim>
    class Source : public dealii::Function<dim>
    {
    public:
      /**
       * @brief Constructor
       *
       * @param physical_parameters User defined runtime parameters
       * @param exp_order Maximum expansion order \f$ l_{\rm max} \f$
       */
      Source(const PhysicalParameters &physical_parameters,
             unsigned int              exp_order)
        : dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
        , prm{physical_parameters}
        , lms_indices{PDESystem::create_lms_indices(exp_order)}
      {}


      /**
       * @brief Evaluate the source term
       *
       * Return a vector of values corresponding to the spherical harmonics
       * decomposition of the source term \f$ S_{i(l,m,s)}(t, \mathbf{x}, p)
       * \f$.
       *
       * @param p Point in reduced phase space
       * @param values Return vector \f$ S_{i(l,m,s)}(t, \mathbf{x}, p) \f$
       * @see InitialValueFunction::vector_value()
       * @see @dealref{Function::vector_value(),classFunction,ae316ebc05d21989d573024f8a23c49cb}
       */
      void
      vector_value(const dealii::Point<dim> &p,
                   dealii::Vector<double>   &values) const override
      {
        AssertDimension(values.size(), this->n_components);
        static_cast<void>(p); // suppress compiler warning

        for (unsigned int i = 0; i < values.size(); ++i)
          {
            /** [Source] */
            // !!!EDIT HERE!!!
            values[i] = 0.;
            /** [Source] */
          }
      }



    private:
      /** User defined runtime parameters */
      const PhysicalParameters prm;
      /**
       * Map between system index \f$ i \f$ and spherical harmonic indices
       * \f$ (l,m,n) \f$
       */
      const std::vector<std::array<unsigned int, 3>> lms_indices;
    };



    /**
     * @brief Magnetic field
     *
     * @tparam dim Dimension of the reduced phase space \f$ (\mathbf{x}, p) \f$
     */
    template <unsigned int dim>
    class MagneticField : public dealii::Function<dim>
    {
    public:
      /**
       * @brief Constructor
       *
       * @param physical_parameters User defined runtime parameters
       */
      MagneticField(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(3)
        , prm{physical_parameters}
      {}


      /**
       * @brief Evaluate the magnetic field at a point `p`
       *
       * Return a vector of values corresponding to the magnetic field
       * \f$ \mathbf{B}(t, \mathbf{x}) \f$ at point in space.
       *
       * Be aware that `point` is a point in reduced phase space \f$
       * (\mathbf{x}, p) \f$ and **not** in configuration space \f$ (\mathbf{x})
       * \f$. Therefore, the first `dim_cs` components correspond to the spatial
       * coordinates \f$ \mathbf{x} \f$, while the last component corresponds to
       * the momentum coordinate \f$ p \f$ and should not be used here (if the
       * momentum term is activated).
       *
       * ```cpp
       * const double x = point[0];
       * const double y = point[1]; // only if dim_cs > 1
       * const double z = point[2]; // only if dim_cs > 2
       * ```
       *
       * The magnetic field has always three components, independent of the
       * dimension of the configuration space. This allows to describe a
       * magnetic field that is points out of the simulation plane.
       *
       * ```cpp
       * magnetic_field[0] = B_x;
       * magnetic_field[1] = B_y;
       * magnetic_field[2] = B_z;
       * ```
       *
       * @param p Point in reduced phase space
       * @param magnetic_field Return vector \f$ \mathbf{B}(t, \mathbf{x}) \f$
       * @see @dealref{Function::vector_value(),classFunction,ae316ebc05d21989d573024f8a23c49cb}
       */
      void
      vector_value(const dealii::Point<dim> &p,
                   dealii::Vector<double>   &magnetic_field) const override
      {
        AssertDimension(magnetic_field.size(), this->n_components);
        static_cast<void>(p); // suppress compiler warning

        /** [Magnetic field] */
        // !!!EDIT HERE!!!
        magnetic_field[0] = 0.;
        magnetic_field[1] = 0.;
        magnetic_field[2] = 0.;
        /** [Magnetic field] */
      }



    private:
      /** User defined runtime parameters */
      const PhysicalParameters prm;
    };



    /**
     * @brief Background velocity field
     *
     * @tparam dim Dimension of the reduced phase space \f$ (\mathbf{x}, p) \f$
     */
    template <unsigned int dim>
    class BackgroundVelocityField : public dealii::Function<dim>
    {
    public:
      /**
       * @brief Constructor
       *
       * @param physical_parameters User defined runtime parameters
       */
      BackgroundVelocityField(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(3)
        , prm{physical_parameters}
      {}



      /**
       * @brief Evaluate the velocity field at a point `p`
       *
       * Return a vector of values corresponding to the velocity field
       * \f$ \mathbf{u}(t, \mathbf{x}) \f$ at point in space.
       *
       * @param p Point in reduced phase space
       * @param velocity Return vector \f$ \mathbf{u}(t, \mathbf{x}) \f$
       * @see MagneticField::vector_value()
       * @see @dealref{Function::vector_value(),classFunction,ae316ebc05d21989d573024f8a23c49cb}
       */
      void
      vector_value(const dealii::Point<dim> &p,
                   dealii::Vector<double>   &velocity) const override
      {
        AssertDimension(velocity.size(), this->n_components);
        static_cast<void>(p); // suppress compiler warning

        /** [Background velocity value] */
        // !!!EDIT HERE!!!
        velocity[0] = 0.;
        velocity[1] = 0.;
        velocity[2] = 0.;
        /** [Background velocity value] */
      }


      /**
       * @brief Evaluate the divergence of the velocity field
       *
       * Return a vector of values corresponding to the divergence of the
       * velocity field \f$ \nabla \cdot \mathbf{u}(t, \mathbf{x}) \f$ at all
       * `points` in space.
       *
       * @param points List of points in reduced phase space
       * @param divergence List of divergence values
       * @see ScatteringFrequency::value_list()
       * @see MagneticField::vector_value()
       * @see @dealref{Function::value_list(),classFunction,abe86ee7f7f12cf4041d1e714c0fb42f3}
       */
      void
      divergence_list(const std::vector<dealii::Point<dim>> &points,
                      std::vector<double>                   &divergence)
      {
        AssertDimension(divergence.size(), points.size());
        static_cast<void>(points); // suppress compiler warning

        for (unsigned int i = 0; i < divergence.size(); ++i)
          {
            /** [Background velocity divergence] */
            // !!!EDIT HERE!!!
            divergence[i] = 0.;
            /** [Background velocity divergence] */
          }
      }



      /**
       * @brief Evaluate the material derivative of the velocity field
       *
       * Return a vector of values corresponding to the material derivative of
       * the velocity field \f$ \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t}
       * \f$ at all `points` in space.
       *
       * @param points List of points in reduced phase space
       * @param material_derivatives List of material derivative values
       * @see ScatteringFrequency::value_list()
       * @see MagneticField::vector_value()
       * @see @dealref{Function::value_list(),classFunction,abe86ee7f7f12cf4041d1e714c0fb42f3}
       *
       */
      void
      material_derivative_list(
        const std::vector<dealii::Point<dim>> &points,
        std::vector<dealii::Vector<double>>   &material_derivatives)
      {
        AssertDimension(material_derivatives.size(), points.size());
        AssertDimension(material_derivatives[0].size(), this->n_components);

        for (unsigned int i = 0; i < points.size(); ++i)
          {
            /** [Background velocity material derivative] */
            // !!!EDIT HERE!!!
            material_derivatives[i][0] = 0.; // D/Dt u_x
            material_derivatives[i][1] = 0.; // D/Dt u_y
            material_derivatives[i][2] = 0.; // D/Dt u_z
            /** [Background velocity material derivative] */
          }
      }



      /**
       * @brief Evaluate the Jacobian matrix of the velocity field
       *
       * Return a vector of matrices corresponding to the Jacobian matrix of
       * the velocity field \f$ \nabla \mathbf{u}(t, \mathbf{x}) \f$ at all
       * `points` in space.
       *
       * @todo Add more documentation, maybe change to vector of FullMatrix?
       *
       * @param points List of points in reduced phase space
       * @param jacobians List of Jacobian matrices
       * @see ScatteringFrequency::value_list()
       * @see MagneticField::vector_value()
       * @see @dealref{Function::value_list(),classFunction,abe86ee7f7f12cf4041d1e714c0fb42f3}
       */
      void
      jacobian_list(
        const std::vector<dealii::Point<dim>>            &points,
        std::vector<std::vector<dealii::Vector<double>>> &jacobians) const
      {
        AssertDimension(jacobians.size(), points.size());
        AssertDimension(jacobians[0].size(), this->n_components);
        AssertDimension(jacobians[0][0].size(), this->n_components);

        for (unsigned int i = 0; i < points.size(); ++i)
          {
            /** [Background velocity Jacobian] */
            // !!!EDIT HERE!!!
            jacobians[i][0][0] = 0.; // \partial u_x / \partial x
            jacobians[i][0][1] = 0.; // \partial u_x / \partial y
            jacobians[i][0][2] = 0.; // \partial u_x / \partial z

            jacobians[i][1][0] = 0.; // \partial u_y / \partial x
            jacobians[i][1][1] = 0.; // \partial u_y / \partial y
            jacobians[i][1][2] = 0.; // \partial u_y / \partial z

            jacobians[i][2][0] = 0.; // \partial u_z / \partial x
            jacobians[i][2][1] = 0.; // \partial u_z / \partial y
            jacobians[i][2][2] = 0.; // \partial u_z / \partial z
            /** [Background velocity Jacobian] */
          }
      }



    private:
      /** User defined runtime parameters */
      const PhysicalParameters prm;
    };
  } // namespace VFP
} // namespace sapphirepp
#endif
