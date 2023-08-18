#include "euler-operator.h"


template <int dim, int degree, int n_points_1d>
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::EulerOperator(
  TimerOutput &timer)
  : timer(timer)
{}



template <int dim, int degree, int n_points_1d>
void
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::reinit(
  const Mapping<dim>    &mapping,
  const DoFHandler<dim> &dof_handler)
{
  const std::vector<const DoFHandler<dim> *> dof_handlers = {&dof_handler};
  const AffineConstraints<double>            dummy;
  const std::vector<const AffineConstraints<double> *> constraints = {&dummy};
  const std::vector<Quadrature<1>> quadratures = {QGauss<1>(n_q_points_1d),
                                                  QGauss<1>(fe_degree + 1)};

  typename MatrixFree<dim, Number>::AdditionalData additional_data;
  additional_data.mapping_update_flags =
    (update_gradients | update_JxW_values | update_quadrature_points |
     update_values);
  additional_data.mapping_update_flags_inner_faces =
    (update_JxW_values | update_quadrature_points | update_normal_vectors |
     update_values);
  additional_data.mapping_update_flags_boundary_faces =
    (update_JxW_values | update_quadrature_points | update_normal_vectors |
     update_values);
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, Number>::AdditionalData::none;

  data.reinit(mapping, dof_handlers, constraints, quadratures, additional_data);
}



template <int dim, int degree, int n_points_1d>
void
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::initialize_vector(
  LinearAlgebra::distributed::Vector<Number> &vector) const
{
  data.initialize_dof_vector(vector);
}



template <int dim, int degree, int n_points_1d>
void
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::set_inflow_boundary(
  const types::boundary_id       boundary_id,
  std::unique_ptr<Function<dim>> inflow_function)
{
  AssertThrow(subsonic_outflow_boundaries.find(boundary_id) ==
                  subsonic_outflow_boundaries.end() &&
                wall_boundaries.find(boundary_id) == wall_boundaries.end(),
              ExcMessage("You already set the boundary with id " +
                         std::to_string(static_cast<int>(boundary_id)) +
                         " to another type of boundary before now setting " +
                         "it as inflow"));
  AssertThrow(inflow_function->n_components == dim + 2,
              ExcMessage("Expected function with dim+2 components"));

  inflow_boundaries[boundary_id] = std::move(inflow_function);
}


template <int dim, int degree, int n_points_1d>
void
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::
  set_subsonic_outflow_boundary(const types::boundary_id       boundary_id,
                                std::unique_ptr<Function<dim>> outflow_function)
{
  AssertThrow(inflow_boundaries.find(boundary_id) == inflow_boundaries.end() &&
                wall_boundaries.find(boundary_id) == wall_boundaries.end(),
              ExcMessage("You already set the boundary with id " +
                         std::to_string(static_cast<int>(boundary_id)) +
                         " to another type of boundary before now setting " +
                         "it as subsonic outflow"));
  AssertThrow(outflow_function->n_components == dim + 2,
              ExcMessage("Expected function with dim+2 components"));

  subsonic_outflow_boundaries[boundary_id] = std::move(outflow_function);
}


template <int dim, int degree, int n_points_1d>
void
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::set_wall_boundary(
  const types::boundary_id boundary_id)
{
  AssertThrow(inflow_boundaries.find(boundary_id) == inflow_boundaries.end() &&
                subsonic_outflow_boundaries.find(boundary_id) ==
                  subsonic_outflow_boundaries.end(),
              ExcMessage("You already set the boundary with id " +
                         std::to_string(static_cast<int>(boundary_id)) +
                         " to another type of boundary before now setting " +
                         "it as wall boundary"));

  wall_boundaries.insert(boundary_id);
}


template <int dim, int degree, int n_points_1d>
void
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::set_body_force(
  std::unique_ptr<Function<dim>> body_force)
{
  AssertDimension(body_force->n_components, dim);

  this->body_force = std::move(body_force);
}



template <int dim, int degree, int n_points_1d>
void
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::local_apply_cell(
  const MatrixFree<dim, Number> &,
  LinearAlgebra::distributed::Vector<Number>       &dst,
  const LinearAlgebra::distributed::Vector<Number> &src,
  const std::pair<unsigned int, unsigned int>      &cell_range) const
{
  FEEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data);

  Tensor<1, dim, VectorizedArray<Number>> constant_body_force;
  const Functions::ConstantFunction<dim> *constant_function =
    dynamic_cast<Functions::ConstantFunction<dim> *>(body_force.get());

  if (constant_function)
    constant_body_force = evaluate_function<dim, Number, dim>(
      *constant_function, Point<dim, VectorizedArray<Number>>());

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      phi.reinit(cell);
      phi.gather_evaluate(src, EvaluationFlags::values);

      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        {
          const auto w_q = phi.get_value(q);
          phi.submit_gradient(euler_flux<dim>(w_q), q);
          if (body_force.get() != nullptr)
            {
              const Tensor<1, dim, VectorizedArray<Number>> force =
                constant_function ?
                  constant_body_force :
                  evaluate_function<dim, Number, dim>(*body_force,
                                                      phi.quadrature_point(q));

              Tensor<1, dim + 2, VectorizedArray<Number>> forcing;
              for (unsigned int d = 0; d < dim; ++d)
                forcing[d + 1] = w_q[0] * force[d];
              for (unsigned int d = 0; d < dim; ++d)
                forcing[dim + 1] += force[d] * w_q[d + 1];

              phi.submit_value(forcing, q);
            }
        }

      phi.integrate_scatter(((body_force.get() != nullptr) ?
                               EvaluationFlags::values :
                               EvaluationFlags::nothing) |
                              EvaluationFlags::gradients,
                            dst);
    }
}



template <int dim, int degree, int n_points_1d>
void
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::local_apply_face(
  const MatrixFree<dim, Number> &,
  LinearAlgebra::distributed::Vector<Number>       &dst,
  const LinearAlgebra::distributed::Vector<Number> &src,
  const std::pair<unsigned int, unsigned int>      &face_range) const
{
  FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi_m(data, true);
  FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi_p(data,
                                                                    false);

  for (unsigned int face = face_range.first; face < face_range.second; ++face)
    {
      phi_p.reinit(face);
      phi_p.gather_evaluate(src, EvaluationFlags::values);

      phi_m.reinit(face);
      phi_m.gather_evaluate(src, EvaluationFlags::values);

      for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
        {
          const auto numerical_flux =
            euler_numerical_flux<dim>(phi_m.get_value(q),
                                      phi_p.get_value(q),
                                      phi_m.get_normal_vector(q));
          phi_m.submit_value(-numerical_flux, q);
          phi_p.submit_value(numerical_flux, q);
        }

      phi_p.integrate_scatter(EvaluationFlags::values, dst);
      phi_m.integrate_scatter(EvaluationFlags::values, dst);
    }
}



template <int dim, int degree, int n_points_1d>
void
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::
  local_apply_boundary_face(
    const MatrixFree<dim, Number> &,
    LinearAlgebra::distributed::Vector<Number>       &dst,
    const LinearAlgebra::distributed::Vector<Number> &src,
    const std::pair<unsigned int, unsigned int>      &face_range) const
{
  FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data, true);

  for (unsigned int face = face_range.first; face < face_range.second; ++face)
    {
      phi.reinit(face);
      phi.gather_evaluate(src, EvaluationFlags::values);

      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        {
          const auto w_m    = phi.get_value(q);
          const auto normal = phi.get_normal_vector(q);

          auto rho_u_dot_n = w_m[1] * normal[0];
          for (unsigned int d = 1; d < dim; ++d)
            rho_u_dot_n += w_m[1 + d] * normal[d];

          bool at_outflow = false;

          Tensor<1, dim + 2, VectorizedArray<Number>> w_p;
          const auto boundary_id = data.get_boundary_id(face);
          if (wall_boundaries.find(boundary_id) != wall_boundaries.end())
            {
              w_p[0] = w_m[0];
              for (unsigned int d = 0; d < dim; ++d)
                w_p[d + 1] = w_m[d + 1] - 2. * rho_u_dot_n * normal[d];
              w_p[dim + 1] = w_m[dim + 1];
            }
          else if (inflow_boundaries.find(boundary_id) !=
                   inflow_boundaries.end())
            w_p =
              evaluate_function(*inflow_boundaries.find(boundary_id)->second,
                                phi.quadrature_point(q));
          else if (subsonic_outflow_boundaries.find(boundary_id) !=
                   subsonic_outflow_boundaries.end())
            {
              w_p          = w_m;
              w_p[dim + 1] = evaluate_function(
                *subsonic_outflow_boundaries.find(boundary_id)->second,
                phi.quadrature_point(q),
                dim + 1);
              at_outflow = true;
            }
          else
            AssertThrow(false,
                        ExcMessage("Unknown boundary id, did "
                                   "you set a boundary condition for "
                                   "this part of the domain boundary?"));

          auto flux = euler_numerical_flux<dim>(w_m, w_p, normal);

          if (at_outflow)
            for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
              {
                if (rho_u_dot_n[v] < -1e-12)
                  for (unsigned int d = 0; d < dim; ++d)
                    flux[d + 1][v] = 0.;
              }

          phi.submit_value(-flux, q);
        }

      phi.integrate_scatter(EvaluationFlags::values, dst);
    }
}



template <int dim, int degree, int n_points_1d>
void
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::
  local_apply_inverse_mass_matrix(
    const MatrixFree<dim, Number> &,
    LinearAlgebra::distributed::Vector<Number>       &dst,
    const LinearAlgebra::distributed::Vector<Number> &src,
    const std::pair<unsigned int, unsigned int>      &cell_range) const
{
  FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1);
  MatrixFreeOperators::CellwiseInverseMassMatrix<dim, degree, dim + 2, Number>
    inverse(phi);

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      phi.reinit(cell);
      phi.read_dof_values(src);

      inverse.apply(phi.begin_dof_values(), phi.begin_dof_values());

      phi.set_dof_values(dst);
    }
}



template <int dim, int degree, int n_points_1d>
void
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::apply(
  const double                                      current_time,
  const LinearAlgebra::distributed::Vector<Number> &src,
  LinearAlgebra::distributed::Vector<Number>       &dst) const
{
  {
    TimerOutput::Scope t(timer, "apply - integrals");

    for (auto &i : inflow_boundaries)
      i.second->set_time(current_time);
    for (auto &i : subsonic_outflow_boundaries)
      i.second->set_time(current_time);

    data.loop(&EulerOperator::local_apply_cell,
              &EulerOperator::local_apply_face,
              &EulerOperator::local_apply_boundary_face,
              this,
              dst,
              src,
              true,
              MatrixFree<dim, Number>::DataAccessOnFaces::values,
              MatrixFree<dim, Number>::DataAccessOnFaces::values);
  }

  {
    TimerOutput::Scope t(timer, "apply - inverse mass");

    data.cell_loop(&EulerOperator::local_apply_inverse_mass_matrix,
                   this,
                   dst,
                   dst);
  }
}



template <int dim, int degree, int n_points_1d>
void
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::perform_stage(
  const Number                                      current_time,
  const Number                                      factor_solution,
  const Number                                      factor_ai,
  const LinearAlgebra::distributed::Vector<Number> &current_ri,
  LinearAlgebra::distributed::Vector<Number>       &vec_ki,
  LinearAlgebra::distributed::Vector<Number>       &solution,
  LinearAlgebra::distributed::Vector<Number>       &next_ri) const
{
  {
    TimerOutput::Scope t(timer, "rk_stage - integrals L_h");

    for (auto &i : inflow_boundaries)
      i.second->set_time(current_time);
    for (auto &i : subsonic_outflow_boundaries)
      i.second->set_time(current_time);

    data.loop(&EulerOperator::local_apply_cell,
              &EulerOperator::local_apply_face,
              &EulerOperator::local_apply_boundary_face,
              this,
              vec_ki,
              current_ri,
              true,
              MatrixFree<dim, Number>::DataAccessOnFaces::values,
              MatrixFree<dim, Number>::DataAccessOnFaces::values);
  }


  {
    TimerOutput::Scope t(timer, "rk_stage - inv mass + vec upd");
    data.cell_loop(
      &EulerOperator::local_apply_inverse_mass_matrix,
      this,
      next_ri,
      vec_ki,
      std::function<void(const unsigned int, const unsigned int)>(),
      [&](const unsigned int start_range, const unsigned int end_range) {
        const Number ai = factor_ai;
        const Number bi = factor_solution;
        if (ai == Number())
          {
            //! Remove SIMD pragma if not working
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (unsigned int i = start_range; i < end_range; ++i)
              {
                const Number k_i          = next_ri.local_element(i);
                const Number sol_i        = solution.local_element(i);
                solution.local_element(i) = sol_i + bi * k_i;
              }
          }
        else
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (unsigned int i = start_range; i < end_range; ++i)
              {
                const Number k_i          = next_ri.local_element(i);
                const Number sol_i        = solution.local_element(i);
                solution.local_element(i) = sol_i + bi * k_i;
                next_ri.local_element(i)  = sol_i + ai * k_i;
              }
          }
      });
  }
}



template <int dim, int degree, int n_points_1d>
void
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::project(
  const Function<dim>                        &function,
  LinearAlgebra::distributed::Vector<Number> &solution) const
{
  FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1);
  MatrixFreeOperators::CellwiseInverseMassMatrix<dim, degree, dim + 2, Number>
    inverse(phi);
  solution.zero_out_ghost_values();
  for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
    {
      phi.reinit(cell);
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        phi.submit_dof_value(evaluate_function(function,
                                               phi.quadrature_point(q)),
                             q);
      inverse.transform_from_q_points_to_basis(dim + 2,
                                               phi.begin_dof_values(),
                                               phi.begin_dof_values());
      phi.set_dof_values(solution);
    }
}



template <int dim, int degree, int n_points_1d>
std::array<double, 3>
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::compute_errors(
  const Function<dim>                              &function,
  const LinearAlgebra::distributed::Vector<Number> &solution) const
{
  TimerOutput::Scope t(timer, "compute errors");
  double             errors_squared[3] = {};
  FEEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data, 0, 0);

  for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
    {
      phi.reinit(cell);
      phi.gather_evaluate(solution, EvaluationFlags::values);
      VectorizedArray<Number> local_errors_squared[3] = {};
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        {
          const auto error =
            evaluate_function(function, phi.quadrature_point(q)) -
            phi.get_value(q);
          const auto JxW = phi.JxW(q);

          local_errors_squared[0] += error[0] * error[0] * JxW;
          for (unsigned int d = 0; d < dim; ++d)
            local_errors_squared[1] += (error[d + 1] * error[d + 1]) * JxW;
          local_errors_squared[2] += (error[dim + 1] * error[dim + 1]) * JxW;
        }
      for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
           ++v)
        for (unsigned int d = 0; d < 3; ++d)
          errors_squared[d] += local_errors_squared[d][v];
    }

  Utilities::MPI::sum(errors_squared, MPI_COMM_WORLD, errors_squared);

  std::array<double, 3> errors;
  for (unsigned int d = 0; d < 3; ++d)
    errors[d] = std::sqrt(errors_squared[d]);

  return errors;
}



template <int dim, int degree, int n_points_1d>
double
Sapphire::Hydro::EulerOperator<dim, degree, n_points_1d>::
  compute_cell_transport_speed(
    const LinearAlgebra::distributed::Vector<Number> &solution) const
{
  TimerOutput::Scope t(timer, "compute transport speed");
  Number             max_transport = 0;
  FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1);

  for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
    {
      phi.reinit(cell);
      phi.gather_evaluate(solution, EvaluationFlags::values);
      VectorizedArray<Number> local_max = 0.;
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        {
          const auto solution = phi.get_value(q);
          const auto velocity = euler_velocity<dim>(solution);
          const auto pressure = euler_pressure<dim>(solution);

          const auto inverse_jacobian = phi.inverse_jacobian(q);
          const auto convective_speed = inverse_jacobian * velocity;
          VectorizedArray<Number> convective_limit = 0.;
          for (unsigned int d = 0; d < dim; ++d)
            convective_limit =
              std::max(convective_limit, std::abs(convective_speed[d]));

          const auto speed_of_sound =
            std::sqrt(gamma * pressure * (1. / solution[0]));

          Tensor<1, dim, VectorizedArray<Number>> eigenvector;
          for (unsigned int d = 0; d < dim; ++d)
            eigenvector[d] = 1.;
          for (unsigned int i = 0; i < 5; ++i)
            {
              eigenvector =
                transpose(inverse_jacobian) * (inverse_jacobian * eigenvector);
              VectorizedArray<Number> eigenvector_norm = 0.;
              for (unsigned int d = 0; d < dim; ++d)
                eigenvector_norm =
                  std::max(eigenvector_norm, std::abs(eigenvector[d]));
              eigenvector /= eigenvector_norm;
            }
          const auto jac_times_ev   = inverse_jacobian * eigenvector;
          const auto max_eigenvalue = std::sqrt((jac_times_ev * jac_times_ev) /
                                                (eigenvector * eigenvector));
          local_max =
            std::max(local_max,
                     max_eigenvalue * speed_of_sound + convective_limit);
        }

      for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
           ++v)
        for (unsigned int d = 0; d < 3; ++d)
          max_transport = std::max(max_transport, local_max[v]);
    }

  max_transport = Utilities::MPI::max(max_transport, MPI_COMM_WORLD);

  return max_transport;
}


// explicit instantiation
// TODO_HD: Make general
template class Sapphire::Hydro::EulerOperator<Sapphire::Hydro::dimension,
                                              Sapphire::Hydro::fe_degree,
                                              Sapphire::Hydro::n_q_points_1d>;