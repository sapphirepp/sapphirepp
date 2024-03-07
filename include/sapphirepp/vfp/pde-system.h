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

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <array>
#include <ostream>
#include <vector>

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



      /** @{ */
      const std::vector<dealii::LAPACKFullMatrix<double>> &
      get_advection_matrices() const;
      const std::vector<dealii::LAPACKFullMatrix<double>> &
      get_generator_rotation_matrices() const;
      const dealii::Vector<double> &
      get_collision_matrix() const;
      const std::vector<dealii::LAPACKFullMatrix<double>> &
      get_adv_mat_products() const;
      const std::vector<dealii::LAPACKFullMatrix<double>> &
      get_adv_cross_gen() const;
      const std::vector<dealii::LAPACKFullMatrix<double>> &
      get_t_matrices() const;
      /** @} */



      /** @{ */
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

      void
      print_advection_matrices(std::ostream &os) const;
      void
      print_generator_rotation_matrices(std::ostream &os) const;
      void
      print_collision_matrix(std::ostream &os) const;
      void
      print_adv_mat_products(std::ostream &os) const;
      void
      print_adv_cross_gen(std::ostream &os) const;
      void
      print_t_matrices(std::ostream &os) const;
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
