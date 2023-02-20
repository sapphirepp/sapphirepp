#include "upwind-flux.h"

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

template <int dim>
VFPEquation::UpwindFlux<dim>::UpwindFlux(const PDESystem &system)
    : pde_system{system},
      matrix_size(pde_system.system_size()),
      advection_matrices(3),
      adv_mat_products(6),
      matrix_sum(matrix_size * matrix_size),
      eigenvalues(matrix_size),
      eigenvectors(matrix_size * matrix_size),
      positive_eigenvalues(matrix_size),
      negative_eigenvalues(matrix_size),
      eigenvectors_advection_matrices(3),
      eigenvalues_advection_matrices(matrix_size),
      lwork{-1},
      liwork{-1} {
  // NOTE: Since we very often call compute_matrix_sum and the matrixes classes
  // of dealii do not allow unchecked access to there raw data, we create copies
  // of the matrices in the hope that this additional memory consumption is made
  // up for by perfomance gains.
  const std::vector<dealii::LAPACKFullMatrix<double>> &adv_mat =
      pde_system.get_advection_matrices();
  for (unsigned int k = 0; k < 3; ++k) {
    for (unsigned int i = 0; i < matrix_size; ++i)
      for (unsigned int j = 0; j < matrix_size; ++j)
        advection_matrices[k][i * matrix_size + j] = adv_mat[k](i, j);
  }
  const std::vector<dealii::LAPACKFullMatrix<double>> &adv_mat_prod =
      pde_system.get_adv_mat_products();
  for (unsigned int k = 0; k < 6; ++k) {
    for (unsigned int i = 0; i < matrix_size; ++i)
      for (unsigned int j = 0; j < matrix_size; ++j)
        adv_mat_products[k][i * matrix_size + j] = adv_mat_prod[k](i, j);
  }
}

template <int dim>
void VFPEquation::UpwindFlux<dim>::compute_matrix_sum(
    const double n_p, const double p, const double gamma,
    const dealii::Vector<double> &material_derivative,
    const std::vector<dealii::Vector<double>> &jacobian) {
  for (unsigned int j = 0; j < matrix_size; ++j)
    for (unsigned int i = 0; i < matrix_size; ++i)
      // NOTE: When computing the eigenvalues of matrix_sum, I call a Fortran
      // routine which requires the matrix elements to be stored in column-major
      // order. But since the resulting matrix is symmteric, there is no
      // difference between column-major order and row-major order. I stick with
      // the C convention and use row-major order.
      matrix_sum[i * matrix_size + j] =
          -n_p *
          (gamma * (material_derivative[0] *
                        advection_matrices[0][i * matrix_size + j] +
                    material_derivative[1] *
                        advection_matrices[1][i * matrix_size + j] +
                    material_derivative[2] *
                        advection_matrices[2][i * matrix_size + j]) +
           p * (jacobian[0][0] * adv_mat_products[0][i * matrix_size + j] +
                jacobian[1][1] * adv_mat_products[3][i * matrix_size + j] +
                jacobian[2][2] * adv_mat_products[5][i * matrix_size + j] +
                (jacobian[0][1] + jacobian[1][0]) *
                    adv_mat_products[1][i * matrix_size + j] +
                (jacobian[0][2] + jacobian[2][0]) *
                    adv_mat_products[2][i * matrix_size + j] +
                (jacobian[1][2] + jacobian[2][1]) *
                    adv_mat_products[4][i * matrix_size + j]));
  // matrix_sum[i * matrix_size + j] =
  //     -n_p *
  //     (gamma * (material_derivative[0] * advection_matrices[0](i, j) +
  //               material_derivative[1] * advection_matrices[1](i, j) +
  //               material_derivative[2] * advection_matrices[2](i, j)) +
  //      p * (jacobian[0][0] * adv_mat_products[0](i, j) +
  //           jacobian[1][1] * adv_mat_products[3](i, j) +
  //           jacobian[2][2] * adv_mat_products[5](i, j) +
  //           (jacobian[0][1] + jacobian[1][0]) * adv_mat_products[1](i, j) +
  //           (jacobian[0][2] + jacobian[2][0]) * adv_mat_products[2](i, j) +
  //           (jacobian[1][2] + jacobian[2][1]) * adv_mat_products[4](i, j)));
}
