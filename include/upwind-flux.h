#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/vector.h>

#include <vector>

#include "particle-functions.h"
#include "pde-system.h"
#include "physical-setup.h"

namespace VFPEquation {
template <int dim>
class UpwindFlux {
  UpwindFlux(const PDESystem& system);
  void set_time();
  void compute_upwind_fluxes(
      const std::vector<dealii::Point<dim>>& q_points,
      const std::vector<dealii::Tensor<1, dim>>& normals,
      std::vector<dealii::FullMatrix<double>>& positive_flux_matrices,
      std::vector<dealii::FullMatrix<double>>& negative_flux_matrices);

 private:
  void compute_eigenvalues_and_eigenvectors(const char range, const double vl,
                                            const double vu,
                                            const double abstol);
  void prepare_upwind_fluxes();
  void separate_positive_and_negative_eigenvalues();
  void triple_product(dealii::FullMatrix<double>& flux_matrix);
  void compute_flux_in_space_directions();
  void compute_matrix_sum(const double n_p, const double p, const double gamma,
                          const dealii::Vector<double>& material_derivative,
                          const std::vector<dealii::Vector<double>>& jacobian);
  void compute_flux_in_p_direction();
  // The computation of the upwind flux in the p direction requires to the
  // compute eigenvalues and eigenvectors for every face with n_p = 1 (or -1).
  // Lapack needs correctly sized work arrays to do these computations. It can
  // determine the sizes it needs with a specific call to the used eigenvalue
  // routine. This is done in constructor of this class. The work arrays will be
  // allocated only once and then reused every time eigenvalues and eigenvectors
  // are computed.
  dealii::types::blas_int info;  // Lapack's error code
  // working array containting doubles
  std::vector<double> work;
  // length of the working array -> to be  determined by the first call of syevr
  dealii::types::blas_int lwork;
  // working array containing integers
  std::vector<dealii::types::blas_int> iwork;
  dealii::types::blas_int liwork;

  std::vector<std::vector<double>> advection_matrices;
  std::vector<std::vector<double>> adv_mat_products;

  // Will be overwritten in flux computations
  std::vector<double> matrix_sum;
  std::vector<double> eigenvalues;
  std::vector<double> eigenvectors;
  std::vector<double> positive_eigenvalues;
  std::vector<double> negative_eigenvalues;

  // Are used in the flux computations in space directions
  std::vector<dealii::LAPACKFullMatrix<double>> eigenvectors_advection_matrices;
  std::vector<double> eigenvalues_advection_matrices;

  // Physical input
  BackgroundVelocityField<dim> background_velocity_field;
  ParticleGamma<dim> particle_gammas;

  // PDE System matrices
  const unsigned int matrix_size;
  const PDESystem& pde_system;
};
}  // namespace VFPEquation
