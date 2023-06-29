#ifndef VFPEQUATION_PDESYSTEM_H
#define VFPEQUATION_PDESYSTEM_H

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <ostream>
#include <vector>

namespace Sapphire {
class PDESystem {
 public:
  PDESystem(int l);
  // get functions
  // matrices
  const std::vector<dealii::LAPACKFullMatrix<double>>& get_advection_matrices()
      const;
  const std::vector<dealii::LAPACKFullMatrix<double>>& get_generator_matrices()
      const;
  const dealii::Vector<double>& get_collision_matrix() const;
  const std::vector<dealii::LAPACKFullMatrix<double>>& get_adv_mat_products()
      const;
  const std::vector<dealii::LAPACKFullMatrix<double>>& get_adv_cross_gen()
      const;
  const std::vector<dealii::LAPACKFullMatrix<double>>& get_t_matrices() const;
  // lms indices
  const std::vector<std::array<unsigned int, 3>>& get_lms_indices() const;

  // returns the size of the system
  unsigned int system_size() const;

  // print functions
  // matrices
  void print_advection_matrices(std::ostream& os) const;
  void print_generator_matrices(std::ostream& os) const;
  void print_collision_matrix(std::ostream& os) const;
  void print_adv_mat_products(std::ostream& os) const;
  void print_adv_cross_gen(std::ostream& os) const;
  void print_t_matrices(std::ostream& os) const;
  void print_pde_system(std::ostream& os) const;
  // lms indices
  template <typename StreamType>
  void print_index_map(StreamType& os) const;

 private:
  void create_advection_matrices();
  void create_generator_rotation_matrices();
  void create_collision_matrix();
  void compute_adv_mat_products();
  void compute_adv_cross_generators();
  void compute_t_matrices();
  void shrink_matrices();

  // PDE System data
  int expansion_order;     // Order of the spherical harmonic expansion
  unsigned int system_sz;  // size of system , i.e. of the quadaratic
                           // matrices
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
  // Map between i and l,m,s (implemented in constructor)
  std::vector<std::array<unsigned int, 3>> lms_indices;
};
}  // namespace Sapphire
#endif