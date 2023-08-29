#include "pde-system.h"

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/lac/lapack_full_matrix.h>

Sapphire::VFP::PDESystem::PDESystem(int l)
  : expansion_order{l}
  , system_sz{static_cast<unsigned int>((l + 1) * (l + 1))}
  , advection_matrices(3)
  , generator_rotation_matrices(3)
  , adv_mat_products(6)
  , adv_x_gen_matrices(3)
  , t_matrices(9)
  , lms_indices((l + 1) * (l + 1))
{
  create_advection_matrices();
  create_generator_rotation_matrices();
  create_collision_matrix();
  compute_adv_mat_products();
  compute_adv_cross_generators();
  compute_t_matrices();
  shrink_matrices();

  // Create the map between i and the the lms indices
  for (int s = 0, idx = 0; s <= 1; ++s)
    {
      for (int l = 0; l <= expansion_order; ++l)
        {
          for (int m = l; m >= s; --m)
            {
              idx                 = l * (l + 1) - (s ? -1 : 1) * m;
              lms_indices[idx][0] = l;
              lms_indices[idx][1] = m;
              lms_indices[idx][2] = s;
            }
        }
    }
}

const std::vector<dealii::LAPACKFullMatrix<double>> &
Sapphire::VFP::PDESystem::get_advection_matrices() const
{
  return advection_matrices;
}
const std::vector<dealii::LAPACKFullMatrix<double>> &
Sapphire::VFP::PDESystem::get_generator_matrices() const
{
  return generator_rotation_matrices;
}
const dealii::Vector<double> &
Sapphire::VFP::PDESystem::get_collision_matrix() const
{
  return collision_matrix;
}
const std::vector<dealii::LAPACKFullMatrix<double>> &
Sapphire::VFP::PDESystem::get_adv_mat_products() const
{
  return adv_mat_products;
}
const std::vector<dealii::LAPACKFullMatrix<double>> &
Sapphire::VFP::PDESystem::get_adv_cross_gen() const
{
  return adv_x_gen_matrices;
}
const std::vector<dealii::LAPACKFullMatrix<double>> &
Sapphire::VFP::PDESystem::get_t_matrices() const
{
  return t_matrices;
}

const std::vector<std::array<unsigned int, 3>> &
Sapphire::VFP::PDESystem::get_lms_indices() const
{
  return lms_indices;
}

unsigned int
Sapphire::VFP::PDESystem::system_size() const
{
  return system_sz;
}

void
Sapphire::VFP::PDESystem::print_advection_matrices(std::ostream &os) const
{
  char subscript = 'x';
  for (const auto &advection_matrix : advection_matrices)
    {
      os << "A_" << subscript << ": " << std::endl;
      advection_matrix.print_formatted(os);
      subscript++;
    }
}

void
Sapphire::VFP::PDESystem::print_generator_matrices(std::ostream &os) const
{
  char subscript = 'x';
  for (const auto &generator_matrix : generator_rotation_matrices)
    {
      os << "Omega_" << subscript << ": " << std::endl;
      generator_matrix.print_formatted(os);
      subscript++;
    }
}

void
Sapphire::VFP::PDESystem::print_collision_matrix(std::ostream &os) const
{
  os << "Collision matrix: \n";
  collision_matrix.print(os);
}

void
Sapphire::VFP::PDESystem::print_adv_mat_products(std::ostream &os) const
{
  char subscript_1 = 'x';
  char subscript_2 = 'x';
  for (unsigned int i = 0; i < 3; ++i)
    {
      for (unsigned int j = i; j < 3; ++j)
        {
          os << "A_" << subscript_1 << subscript_2 << ": " << std::endl;
          adv_mat_products[3 * i - i * (i + 1) / 2 + j].print_formatted(os);
          subscript_2++;
        }
      subscript_1++;
      subscript_2 = subscript_1;
    }
}

void
Sapphire::VFP::PDESystem::print_adv_cross_gen(std::ostream &os) const
{
  char subscript = 'x';
  for (const auto &adv_x_gen_mat : adv_x_gen_matrices)
    {
      os << "(A x Omega)_" << subscript << ": " << std::endl;
      adv_x_gen_mat.print_formatted(std::cout);
      subscript++;
    }
}

void
Sapphire::VFP::PDESystem::print_t_matrices(std::ostream &os) const
{
  char subscript_1 = 'x';
  char subscript_2 = 'x';
  for (unsigned int i = 0; i < 3; ++i)
    {
      for (unsigned int j = 0; j < 3; ++j)
        {
          os << "T_" << subscript_1 << subscript_2 << ": " << std::endl;
          t_matrices[3 * i + j].print_formatted(os);
          subscript_2++;
        }
      subscript_2 = 'x';
      subscript_1++;
    }
}

template <typename StreamType>
void
Sapphire::VFP::PDESystem::print_index_map(StreamType &os) const
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

// explicit instantiation
template void
Sapphire::VFP::PDESystem::print_index_map(std::ostream &os) const;
template void
Sapphire::VFP::PDESystem::print_index_map(dealii::ConditionalOStream &os) const;

void
Sapphire::VFP::PDESystem::print_pde_system(std::ostream &os) const
{
  print_advection_matrices(os);
  print_generator_matrices(os);
  print_collision_matrix(os);
  print_adv_mat_products(os);
  print_adv_cross_gen(os);
  print_t_matrices(os);
}

void
Sapphire::VFP::PDESystem::create_advection_matrices()
{
  // The matrix products (e.g. A_x * A_y) only yield the correct system if we
  // compute the matrices for expansion_order + 1 and later shrink them to
  // expansion_order
  unsigned int matrix_size = (expansion_order + 2) * (expansion_order + 2);
  for (auto &advection_matrix : advection_matrices)
    advection_matrix.reinit(matrix_size);

  for (int s = 0; s <= 1; ++s)
    {
      for (int l = 0, i = 0; l <= expansion_order + 1; ++l)
        {
          for (int m = l; m >= s; --m)
            {
              i = l * (l + 1) - (s ? -1 : 1) * m; // (-1)^s
              for (int s_prime = 0; s_prime <= 1; ++s_prime)
                {
                  for (int l_prime = 0, j = 0; l_prime <= expansion_order + 1;
                       ++l_prime)
                    {
                      for (int m_prime = l_prime; m_prime >= s_prime; --m_prime)
                        {
                          j = l_prime * (l_prime + 1) -
                              (s_prime ? -1 : 1) * m_prime;
                          // Ax
                          if (l + 1 == l_prime && m == m_prime && s == s_prime)
                            advection_matrices[0].set(
                              i,
                              j,
                              std::sqrt(((l - m + 1.) * (l + m + 1.)) /
                                        ((2. * l + 3.) * (2 * l + 1.))));
                          if (l - 1 == l_prime && m == m_prime && s == s_prime)
                            advection_matrices[0].set(
                              i,
                              j,
                              std::sqrt(((l - m) * (l + m)) /
                                        ((2. * l + 1.) * (2. * l - 1.))));
                          // Ay
                          if ((l + 1) == l_prime && (m + 1) == m_prime &&
                              s == s_prime)
                            advection_matrices[1].set(
                              i,
                              j,
                              -0.5 *
                                std::sqrt(((l + m + 1.) * (l + m + 2.)) /
                                          ((2. * l + 3.) * (2. * l + 1.))));
                          if ((l - 1) == l_prime && m + 1 == m_prime &&
                              s == s_prime)
                            advection_matrices[1].set(
                              i,
                              j,
                              0.5 * std::sqrt(((l - m - 1.) * (l - m)) /
                                              ((2. * l + 1.) * (2. * l - 1.))));
                          if ((l + 1) == l_prime && (m - 1) == m_prime &&
                              s == s_prime)
                            advection_matrices[1].set(
                              i,
                              j,
                              0.5 * std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                                              ((2. * l + 3.) * (2 * l + 1.))));
                          if ((l - 1) == l_prime && (m - 1) == m_prime &&
                              s == s_prime)
                            advection_matrices[1].set(
                              i,
                              j,
                              -0.5 *
                                std::sqrt(((l + m - 1.) * (l + m)) /
                                          ((2. * l + 1.) * (2. * l - 1.))));
                          // Az
                          // l + 1 = l_prime and s = 1 and s_prime = 0
                          if (l + 1 == l_prime && m + 1 == m_prime && s == 1 &&
                              s_prime == 0)
                            advection_matrices[2].set(
                              i,
                              j,
                              0.5 * (m_prime == 0 ? 1 / std::sqrt(2) : 1.) *
                                std::sqrt(((l + m + 1.) * (l + m + 2.)) /
                                          ((2 * l + 3.) * (2 * l + 1.))));
                          if (l + 1 == l_prime && m - 1 == m_prime && s == 1 &&
                              s_prime == 0)
                            advection_matrices[2].set(
                              i,
                              j,
                              0.5 * (m_prime == 0 ? 1 / std::sqrt(2) : 1.) *
                                std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                                          ((2 * l + 3.) * (2 * l + 1.))));
                          // NOTE: The correction are now directly included (->
                          // new way to compute the real matrices)
                          if (l + 1 == l_prime && m == 0 && m_prime == 1 &&
                              s == 1 && s_prime == 0)
                            advection_matrices[2](i, j) -=
                              0.5 * std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                                              ((2 * l + 3.) * (2 * l + 1.)));
                          if (l + 1 == l_prime && m == 1 && m_prime == 0 &&
                              s == 1 && s_prime == 0)
                            advection_matrices[2](i, j) +=
                              0.5 * 1. / std::sqrt(2) *
                              std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                                        ((2 * l + 3.) * (2 * l + 1.)));

                          // l+ 1 = l_prime and s = 0 and s_prime = 1
                          if (l + 1 == l_prime && m + 1 == m_prime && s == 0 &&
                              s_prime == 1)
                            advection_matrices[2].set(
                              i,
                              j,
                              -0.5 * (m == 0 ? 1 / std::sqrt(2) : 1.) *
                                std::sqrt(((l + m + 1.) * (l + m + 2.)) /
                                          ((2 * l + 3.) * (2 * l + 1.))));
                          if (l + 1 == l_prime && m - 1 == m_prime && s == 0 &&
                              s_prime == 1)
                            advection_matrices[2].set(
                              i,
                              j,
                              -0.5 * (m == 0 ? 1 / std::sqrt(2) : 1.) *
                                std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                                          ((2 * l + 3.) * (2 * l + 1.))));
                          // Corrections
                          if (l + 1 == l_prime && m == 0 && m_prime == 1 &&
                              s == 0 && s_prime == 1)
                            advection_matrices[2](i, j) -=
                              0.5 * 1. / std::sqrt(2) *
                              std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                                        ((2 * l + 3.) * (2 * l + 1.)));
                          if (l + 1 == l_prime && m == 1 && m_prime == 0 &&
                              s == 0 && s_prime == 1)
                            advection_matrices[2](i, j) +=
                              0.5 * std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                                              ((2 * l + 3.) * (2 * l + 1.)));

                          // l - 1 = l_prime and s = 1 and s_prime = 0
                          if (l - 1 == l_prime && m + 1 == m_prime && s == 1 &&
                              s_prime == 0)
                            advection_matrices[2].set(
                              i,
                              j,
                              -0.5 * (m_prime == 0 ? 1 / std::sqrt(2) : 1.) *
                                std::sqrt(((l - m - 1.) * (l - m)) /
                                          ((2 * l + 1.) * (2 * l - 1.))));
                          if (l - 1 == l_prime && m - 1 == m_prime && s == 1 &&
                              s_prime == 0)
                            advection_matrices[2].set(
                              i,
                              j,
                              -0.5 * (m_prime == 0 ? 1 / std::sqrt(2) : 1.) *
                                std::sqrt(((l + m - 1.) * (l + m)) /
                                          ((2 * l + 1.) * (2 * l - 1.))));
                          // Corrections
                          if (l - 1 == l_prime && m == 0 && m_prime == 1 &&
                              s == 1 && s_prime == 0)
                            advection_matrices[2](i, j) +=
                              0.5 * std::sqrt(((l + m - 1.) * (l + m)) /
                                              ((2 * l + 1.) * (2 * l - 1.)));
                          if (l - 1 == l_prime && m == 1 && m_prime == 0 &&
                              s == 1 && s_prime == 0)
                            advection_matrices[2](i, j) -=
                              0.5 * 1. / std::sqrt(2) *
                              std::sqrt(((l + m - 1.) * (l + m)) /
                                        ((2 * l + 1.) * (2 * l - 1.)));

                          // l - 1 = l_prime and s = 0 and s_prime = 1
                          if (l - 1 == l_prime && m + 1 == m_prime && s == 0 &&
                              s_prime == 1)
                            advection_matrices[2].set(
                              i,
                              j,
                              0.5 * (m == 0 ? 1 / std::sqrt(2) : 1.) *
                                std::sqrt(((l - m - 1.) * (l - m)) /
                                          ((2 * l + 1.) * (2 * l - 1.))));
                          if (l - 1 == l_prime && m - 1 == m_prime && s == 0 &&
                              s_prime == 1)
                            advection_matrices[2].set(
                              i,
                              j,
                              0.5 * (m == 0 ? 1 / std::sqrt(2) : 1.) *
                                std::sqrt(((l + m - 1.) * (l + m)) /
                                          ((2 * l + 1.) * (2 * l - 1.))));
                          // Corrections
                          if (l - 1 == l_prime && m == 0 && m_prime == 1 &&
                              s == 0 && s_prime == 1)
                            advection_matrices[2](i, j) +=
                              0.5 * 1. / std::sqrt(2) *
                              std::sqrt(((l + m - 1.) * (l + m)) /
                                        ((2 * l + 1.) * (2 * l - 1.)));
                          if (l - 1 == l_prime && m == 1 && m_prime == 0 &&
                              s == 0 && s_prime == 1)
                            advection_matrices[2](i, j) -=
                              0.5 * std::sqrt(((l + m - 1.) * (l + m)) /
                                              ((2 * l + 1.) * (2 * l - 1.)));
                        }
                    }
                }
            }
        }
    }
  // TODO: Remove the correction part and derive formulas as you did for Az
  for (int l = 0; l <= expansion_order + 1; ++l)
    {
      // Special cases for A_y (necessary for every value of l)
      // Above the diagonal
      // l + 1 = l_prime, m = 0, s = 0, and m_prime = 1, s_prime = 0
      if (l != expansion_order + 1)
        advection_matrices[1](l * (l + 1), (l + 1) * (l + 2) - 1) =
          std::sqrt(2) *
          advection_matrices[1](l * (l + 1), (l + 2) * (l + 1) - 1);
      // l + 1 = l_prime, m = 1, s = 0, and m_prime = 0, s_prime = 0
      if (l != 0 && l != expansion_order + 1)
        advection_matrices[1](l * (l + 1) - 1, (l + 1) * (l + 2)) =
          std::sqrt(2) *
          advection_matrices[1](l * (l + 1) - 1, (l + 2) * (l + 1));
      // Below the diagonal
      // l - 1 = l_prime, m = 0, s = 0, and m_prime = 1, s_prime = 0
      if (l > 1)
        advection_matrices[1](l * (l + 1), l * (l - 1) - 1) =
          std::sqrt(2) * advection_matrices[1](l * (l + 1), l * (l - 1) - 1);
      // l - 1 = l_prime , m = 1, s = 0 and m_prime = 0, s_prime = 0
      if (l != 0)
        advection_matrices[1](l * (l + 1) - 1, l * (l - 1)) =
          std::sqrt(2) * advection_matrices[1](l * (l + 1) - 1, l * (l - 1));
    }
}

void
Sapphire::VFP::PDESystem::create_generator_rotation_matrices()
{
  unsigned int matrix_size = (expansion_order + 2) * (expansion_order + 2);
  for (auto &generator_matrix : generator_rotation_matrices)
    generator_matrix.reinit(matrix_size);

  for (int s = 0; s <= 1; ++s)
    {
      for (int l = 0, i = 0; l <= expansion_order + 1; ++l)
        {
          for (int m = l; m >= s; --m)
            {
              i = l * (l + 1) - (s ? -1 : 1) * m; // (-1)^s
              for (int s_prime = 0; s_prime <= 1; ++s_prime)
                {
                  for (int l_prime = 0, j = 0; l_prime <= expansion_order + 1;
                       ++l_prime)
                    {
                      for (int m_prime = l_prime; m_prime >= s_prime; --m_prime)
                        {
                          j = l_prime * (l_prime + 1) -
                              (s_prime ? -1 : 1) * m_prime;
                          // Omega_x
                          if (l == l_prime && m == m_prime && s == 0 &&
                              s_prime == 1)
                            {
                              generator_rotation_matrices[0].set(i, j, 1. * m);
                              generator_rotation_matrices[0].set(
                                j,
                                i,
                                -1. * m); // Omega matrices are anti-symmetric
                            }
                          // Omega_y
                          if (l == l_prime && (m + 1) == m_prime && s == 0 &&
                              s_prime == 1)
                            {
                              generator_rotation_matrices[1].set(
                                i, j, 0.5 * std::sqrt((l + m + 1.) * (l - m)));
                              generator_rotation_matrices[1].set(
                                j, i, -0.5 * std::sqrt((l + m + 1.) * (l - m)));
                            }
                          if (l == l_prime && (m - 1) == m_prime && s == 0 &&
                              s_prime == 1)
                            {
                              generator_rotation_matrices[1].set(
                                i, j, 0.5 * std::sqrt((l - m + 1.) * (l + m)));
                              generator_rotation_matrices[1].set(
                                j, i, -0.5 * std::sqrt((l - m + 1.) * (l + m)));
                            }
                          // Omega_z
                          if (l == l_prime && (m + 1) == m_prime &&
                              s == s_prime)
                            {
                              generator_rotation_matrices[2].set(
                                i, j, -0.5 * std::sqrt((l + m + 1.) * (l - m)));
                            }
                          if (l == l_prime && (m - 1) == m_prime &&
                              s == s_prime)
                            {
                              generator_rotation_matrices[2].set(
                                i, j, 0.5 * std::sqrt((l - m + 1.) * (l + m)));
                            }
                        }
                    }
                }
            }
        }
    }

  // Edit entries of the matrices around m=0 and s=0. After the transformation
  // they fall out of the pattern of the matrix elements, namely they differ by
  // a factor of 2^(1/2).
  //
  // NOTE: These corrections were not included in the above loops, because some
  // of them were overwritten. This is an effect of the s loops being outside
  // the l and m loops.
  for (int l = 0; l <= expansion_order + 1; ++l)
    {
      // Special cases for Omega
      if (l > 0)
        {
          // l == l_prime, m = 0, s = 0 and m_prime = 1 and s_prime = 1
          generator_rotation_matrices[1](l * (l + 1), l * (l + 1) + 1) =
            std::sqrt(2) *
            generator_rotation_matrices[1](l * (l + 1), l * (l + 1) + 1);
          // l == l_prime, m = 1, s = 1 and m_prime = 0 and s_prime = 0
          generator_rotation_matrices[1](l * (l + 1) + 1, l * (l + 1)) =
            std::sqrt(2) *
            generator_rotation_matrices[1](l * (l + 1) + 1, l * (l + 1));
          // l == l_prime, m = 0, s = 0 and m_prime = 1 and s_prime = 0
          generator_rotation_matrices[2](l * (l + 1), l * (l + 1) - 1) =
            std::sqrt(2) *
            generator_rotation_matrices[2](l * (l + 1), l * (l + 1) - 1);
          // // l == l_prime, m = 1, s = 0 and m_prime = 0 and s_prime = 0
          generator_rotation_matrices[2](l * (l + 1) - 1, l * (l + 1)) =
            std::sqrt(2) *
            generator_rotation_matrices[2](l * (l + 1) - 1, l * (l + 1));
        }
    }
}

void
Sapphire::VFP::PDESystem::create_collision_matrix()
{
  unsigned int matrix_size = (expansion_order + 1) * (expansion_order + 1);
  collision_matrix.reinit(matrix_size);
  for (int s = 0; s <= 1; ++s)
    {
      for (int l = 0, i = 0; l <= expansion_order; ++l)
        {
          for (int m = l; m >= s; --m)
            {
              i = l * (l + 1) - (s ? -1 : 1) * m; // (-1)^s
              for (int s_prime = 0; s_prime <= 1; ++s_prime)
                {
                  for (int l_prime = 0; l_prime <= expansion_order; ++l_prime)
                    {
                      for (int m_prime = l_prime; m_prime >= s_prime; --m_prime)
                        {
                          // C
                          if (l == l_prime && m == m_prime && s == s_prime)
                            {
                              collision_matrix[i] = 0.5 * l * (l + 1.);
                            }
                        }
                    }
                }
            }
        }
    }
}

void
Sapphire::VFP::PDESystem::compute_adv_mat_products()
{
  unsigned int matrix_size = (expansion_order + 2) * (expansion_order + 2);
  for (auto &adv_mat_product : adv_mat_products)
    adv_mat_product.reinit(matrix_size);
  // row-major order, i.e. A_xx, A_xy, A_xz, Ayy, Ayz, Azz
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = i; j < 3; ++j)
      advection_matrices[i].mmult(adv_mat_products[3 * i - i * (i + 1) / 2 + j],
                                  advection_matrices[j]);
}

void
Sapphire::VFP::PDESystem::compute_adv_cross_generators()
{
  unsigned int matrix_size = (expansion_order + 2) * (expansion_order + 2);

  for (auto &adv_x_gen_mat : adv_x_gen_matrices)
    adv_x_gen_mat.reinit(matrix_size);
  // A cross Omega

  // A temporary matrix is needed to add the two matrix products (the interface
  // of dealii of LAPACKFULLMatrix is not well suited to compute the cross
  // product of matrices)
  dealii::LAPACKFullMatrix<double> temp_matrix(matrix_size);

  // A_y * Omega_z - A_z * Omega_y
  advection_matrices[1].mmult(adv_x_gen_matrices[0],
                              generator_rotation_matrices[2]);
  advection_matrices[2].mmult(temp_matrix, generator_rotation_matrices[1]);
  adv_x_gen_matrices[0].add(-1., temp_matrix);
  // A_z * Omega_x - A_x * Omega_z
  advection_matrices[2].mmult(adv_x_gen_matrices[1],
                              generator_rotation_matrices[0]);
  advection_matrices[0].mmult(temp_matrix, generator_rotation_matrices[2]);
  adv_x_gen_matrices[1].add(-1., temp_matrix);
  // A_x * Omega_y - A_y * Omega_x
  advection_matrices[0].mmult(adv_x_gen_matrices[2],
                              generator_rotation_matrices[1]);
  advection_matrices[1].mmult(temp_matrix, generator_rotation_matrices[0]);
  adv_x_gen_matrices[2].add(-1., temp_matrix);
}

void ::Sapphire::VFP::PDESystem::compute_t_matrices()
{
  unsigned int matrix_size = (expansion_order + 2) * (expansion_order + 2);
  for (auto &t_mat : t_matrices)
    t_mat.reinit(matrix_size);
  // the T matrices are stored in row-major order
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      advection_matrices[j].mmult(t_matrices[3 * i + j], adv_x_gen_matrices[i]);
}

void ::Sapphire::VFP::PDESystem::shrink_matrices()
{
  // Shrink the matrices such that they agree with order of the expansion
  unsigned int num_exp_coefficients =
    (expansion_order + 1) * (expansion_order + 1);

  for (auto &advection_matrix : advection_matrices)
    advection_matrix.grow_or_shrink(num_exp_coefficients);

  for (auto &generator_matrix : generator_rotation_matrices)
    generator_matrix.grow_or_shrink(num_exp_coefficients);

  for (auto &adv_mat_product : adv_mat_products)
    adv_mat_product.grow_or_shrink(num_exp_coefficients);

  for (auto &adv_x_gen_mat : adv_x_gen_matrices)
    adv_x_gen_mat.grow_or_shrink(num_exp_coefficients);

  for (auto &t_mat : t_matrices)
    t_mat.grow_or_shrink(num_exp_coefficients);
}
