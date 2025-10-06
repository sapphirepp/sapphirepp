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
 * @file pde-system.h
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::VFP::PDESystem
 */

#ifndef VFP_PDESYSTEM_H
#define VFP_PDESYSTEM_H

#include <deal.II/base/table.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <array>
#include <ostream>
#include <string>
#include <vector>

#include "vfp-flags.h"

namespace sapphirepp
{
  namespace VFP
  {
    /**
     * @brief Calculate the matrices of the PDE system.
     *
     * This class calculates everything related to the PDE system resulting from
     * the VFP equation, including the mapping between the indices and the
     * matrices resulting from the expansion.
     */
    class PDESystem
    {
    public:
      /** @{ */
      /** Order of the spherical harmonic expansion */
      const unsigned int expansion_order;

      /**
       * @brief Size of the system
       *
       * Size of system, i.e. the number of expansion coefficients and the size
       * the quadratic matrices
       * */
      const unsigned int system_size;

      /**
       * @brief Map between system index \f$ i \f$ and spherical harmonic
       *        indices \f$ (l,m,s) \f$
       *
       * The mapping is given by `lms_indices[i] = {l,m,s}`.
       */
      const std::vector<std::array<unsigned int, 3>> lms_indices;
      /** @} */



      /** Constructor */
      PDESystem(const unsigned int expansion_order);


      /** @{ */
      /**
       * @brief Create a mapping between the system index \f$ i \f$ and the
       *        spherical harmonic indices \f$ (l,m,s) \f$
       *
       * @param system_size Number of expansion coefficients, normally
       *        \f$ (l_{\rm max} + 1)^2 \f$.
       * @return std::vector<std::array<unsigned int, 3>> lms_indices.
       *         Mapping `lms_indices[i] = {l,m,s}`.
       */
      static std::vector<std::array<unsigned int, 3>>
      create_lms_indices(const unsigned int system_size);

      /**
       * @brief Create a list of component names, `f_lms`
       *
       * @param system_size Number of expansion coefficients, normally
       *        \f$ (l_{\rm max} + 1)^2 \f$.
       * @param prefix Prefix for the component names.
       * @return std::vector<std::string> component_names.
       *         List `component_names[i] = "f_lms"`.
       */
      static std::vector<std::string>
      create_component_name_list(const unsigned int system_size,
                                 const std::string &prefix = "f_");
      /** @} */



      /**
       * @brief Compute the coupling tables.
       *
       * @param dim_cs Dimension of the configuration space
       * @param vfp_flags VFP Flags
       * @param cell_integrals_mask Coupling table for the cells
       * @param face_integrals_mask Coupling table for the faces
       */
      void
      compute_coupling_tables(
        const unsigned int                            dim_cs,
        const VFPFlags                                vfp_flags,
        dealii::Table<2, dealii::DoFTools::Coupling> &cell_integrals_mask,
        dealii::Table<2, dealii::DoFTools::Coupling> &face_integrals_mask)
        const;



      /** @{ */
      /**
       * @brief Get the advection matrices object
       *
       * @return const std::vector<dealii::LAPACKFullMatrix<double>>&
       */
      const std::vector<dealii::LAPACKFullMatrix<double>> &
      get_advection_matrices() const;

      /**
       * @brief Get the generator rotation matrices object
       *
       * @return const std::vector<dealii::LAPACKFullMatrix<double>>&
       */
      const std::vector<dealii::LAPACKFullMatrix<double>> &
      get_generator_rotation_matrices() const;

      /**
       * @brief Get the collision matrix object
       *
       * @return const dealii::Vector<double>&
       */
      const dealii::Vector<double> &
      get_collision_matrix() const;

      /**
       * @brief Get the adv mat products object
       *
       * @return const std::vector<dealii::LAPACKFullMatrix<double>>&
       */
      const std::vector<dealii::LAPACKFullMatrix<double>> &
      get_adv_mat_products() const;

      /**
       * @brief Get the adv cross gen object
       *
       * @return const std::vector<dealii::LAPACKFullMatrix<double>>&
       */
      const std::vector<dealii::LAPACKFullMatrix<double>> &
      get_adv_cross_gen() const;

      /**
       * @brief Get the t matrices object
       *
       * @return const std::vector<dealii::LAPACKFullMatrix<double>>&
       */
      const std::vector<dealii::LAPACKFullMatrix<double>> &
      get_t_matrices() const;
      /** @} */



      /** @{ */
      /**
       * @brief Print the lms indices ordering
       *
       * @tparam StreamType Type of the output stream
       * @param os Output stream
       */
      template <typename StreamType>
      void
      print_lms_indices(StreamType &os) const
      {
        os << "Ordering of the lms indices: " << std::endl;
        unsigned int i = 0;
        for (const std::array<unsigned int, 3> &lms : lms_indices)
          {
            os << i << ": " << lms[0] << lms[1] << lms[2] << "\n";
            ++i;
          }
        os << std::endl;
      }

      /**
       * @brief Print the advection matrices
       *
       * @param os Output stream
       */
      void
      print_advection_matrices(std::ostream &os) const;

      /**
       * @brief Print the generator rotation matrices
       *
       * @param os Output stream
       */
      void
      print_generator_rotation_matrices(std::ostream &os) const;

      /**
       * @brief Print the collision matrix
       *
       * @param os Output stream
       */
      void
      print_collision_matrix(std::ostream &os) const;

      /**
       * @brief Print the adv mat products
       *
       * @param os Output stream
       */
      void
      print_adv_mat_products(std::ostream &os) const;

      /**
       * @brief Print the adv cross gen
       *
       * @param os Output stream
       */
      void
      print_adv_cross_gen(std::ostream &os) const;

      /**
       * @brief Print the t matrices
       *
       * @param os Output stream
       */
      void
      print_t_matrices(std::ostream &os) const;

      /**
       * @brief Print the full PDE system
       *
       * @param os Output stream
       */
      void
      print_pde_system(std::ostream &os) const;
      /** @} */



    private:
      /** @{ */
      /** Advection matrices */
      std::vector<dealii::LAPACKFullMatrix<double>> advection_matrices;

      /** Rotation matrices (due to the magnetic field) */
      std::vector<dealii::LAPACKFullMatrix<double>> generator_rotation_matrices;

      /** Collision matrix (essentially a reaction matrix) */
      dealii::Vector<double> collision_matrix;

      /** (magnitude) p advection */
      std::vector<dealii::LAPACKFullMatrix<double>> adv_mat_products;

      /** A cross Omega matrices */
      std::vector<dealii::LAPACKFullMatrix<double>> adv_x_gen_matrices;

      /** T matrices */
      std::vector<dealii::LAPACKFullMatrix<double>> t_matrices;
      /** @} */



      void
      create_advection_matrices();
      void
      create_generator_rotation_matrices();
      void
      create_collision_matrix();
      void
      compute_adv_mat_products();
      void
      compute_adv_cross_generators();
      void
      compute_t_matrices();
      void
      shrink_matrices();
    };
  } // namespace VFP
} // namespace sapphirepp
#endif
