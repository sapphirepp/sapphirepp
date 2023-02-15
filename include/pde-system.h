#ifndef VFPEQUATION_PDESYSTEM_H
#define VFPEQUATION_PDESYSTEM_H

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <vector>
namespace VFPEquation {
class PDESystem {
 public:
  PDESystem(int l);
  std::vector<dealii::LAPACKFullMatrix<double>>& get_advection_matrices();
  std::vector<dealii::LAPACKFullMatrix<double>>& get_rotation_matrices();
  std::vector<dealii::LAPACKFullMatrix<double>>& get_adv_mat_products();
  std::vector<dealii::LAPACKFullMatrix<double>>& get_adv_cross_gen();
  std::vector<dealii::LAPACKFullMatrix<double>>& get_t_matrices();

 private:
  void create_advection_matrices();
  void create_generator_rotation_matrices();
  void create_collision_matrix();
  void compute_adv_mat_products();
  void compute_adv_cross_generators();
  void compute_t_matrices();
  void shrink_matrices();

  // PDE System data
  int expansion_order;  // Order of the spherical harmonic expansion
  // Advection matrices
  std::vector<dealii::LAPACKFullMatrix<double>> advection_matrices;
  // Rotation matrices (due to the magnetic field)
  std::vector<dealii::LAPACKFullMatrix<double>> generator_rotation_matrices;
  // (magnitude) p advection
  std::vector<dealii::LAPACKFullMatrix<double>> adv_mat_products;
  // A cross Omega matrices
  std::vector<dealii::LAPACKFullMatrix<double>> adv_x_gen_matrices;
  // T matrices
  std::vector<dealii::LAPACKFullMatrix<double>> t_matrices;
  // Collision matrix (essentially a reaction matrix)
  dealii::Vector<double> collision_matrix;
};
}  // namespace VFPEquation
#endif
