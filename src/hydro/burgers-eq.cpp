#include "burgers-eq.h"

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "numerics.h"

namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;
    template <int dim>
    class ScratchDataDG
    {
    public:
      // Constructor
      ScratchDataDG(
        const Mapping<dim>        &mapping,
        const FiniteElement<dim>  &fe,
        const Quadrature<dim>     &quadrature,
        const Quadrature<dim - 1> &face_quadrature,
        const UpdateFlags update_flags = update_values | update_gradients |
                                         update_quadrature_points |
                                         update_JxW_values,
        const UpdateFlags face_update_flags = update_values |
                                              update_quadrature_points |
                                              update_JxW_values |
                                              update_normal_vectors,
        const UpdateFlags neighbor_face_update_flags = update_values)
        : fe_values(mapping, fe, quadrature, update_flags)
        , fe_face_values(mapping, fe, face_quadrature, face_update_flags)
        , fe_face_values_neighbor(mapping,
                                  fe,
                                  face_quadrature,
                                  neighbor_face_update_flags)
      {}

      // Copy constructor
      ScratchDataDG(const ScratchDataDG &scratch_data)
        : fe_values(scratch_data.fe_values.get_mapping(),
                    scratch_data.fe_values.get_fe(),
                    scratch_data.fe_values.get_quadrature(),
                    scratch_data.fe_values.get_update_flags())
        , fe_face_values(scratch_data.fe_face_values.get_mapping(),
                         scratch_data.fe_face_values.get_fe(),
                         scratch_data.fe_face_values.get_quadrature(),
                         scratch_data.fe_face_values.get_update_flags())
        , fe_face_values_neighbor(
            scratch_data.fe_face_values_neighbor.get_mapping(),
            scratch_data.fe_face_values_neighbor.get_fe(),
            scratch_data.fe_face_values_neighbor.get_quadrature(),
            scratch_data.fe_face_values_neighbor.get_update_flags())
      {}

      FEValues<dim>     fe_values;
      FEFaceValues<dim> fe_face_values;
      FEFaceValues<dim> fe_face_values_neighbor;
    };

    struct CopyDataFaceDG
    {
      Vector<double> cell_vector_1;
      Vector<double> cell_vector_2;

      std::vector<types::global_dof_index> local_dof_indices;
      std::vector<types::global_dof_index> local_dof_indices_neighbor;

      template <typename Iterator>
      void
      reinit(const Iterator &cell,
             const Iterator &neighbor_cell,
             unsigned int    dofs_per_cell)
      {
        cell_vector_1.reinit(dofs_per_cell);
        cell_vector_2.reinit(dofs_per_cell);

        local_dof_indices.resize(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);

        local_dof_indices_neighbor.resize(dofs_per_cell);
        neighbor_cell->get_dof_indices(local_dof_indices_neighbor);
      }
    };

    struct CopyDataDG
    {
      Vector<double>                       cell_vector;
      std::vector<types::global_dof_index> local_dof_indices;
      std::vector<CopyDataFaceDG>          face_data;

      template <typename Iterator>
      void
      reinit(const Iterator &cell, unsigned int dofs_per_cell)
      {
        cell_vector.reinit(dofs_per_cell);

        local_dof_indices.resize(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);
      }
    };

    template <int dim>
    class ScratchDataSlopeLimiter
    {
    public:
      // Constructor
      ScratchDataSlopeLimiter(
        const Mapping<dim>        &mapping,
        const FiniteElement<dim>  &fe,
        const Quadrature<dim>     &quadrature,
        const Quadrature<dim>     &support_quadrature,
        const Quadrature<dim - 1> &face_quadrature,
        const UpdateFlags update_flags = update_values | update_gradients |
                                         update_JxW_values,
        const UpdateFlags neighbor_update_flags = update_values |
                                                  update_JxW_values,
        const UpdateFlags interpolate_update_flags = update_quadrature_points,
        const UpdateFlags face_update_flags = update_values | update_JxW_values)
        : fe_values(mapping, fe, quadrature, update_flags)
        , fe_values_neighbor(mapping, fe, quadrature, neighbor_update_flags)
        , fe_values_interpolate(mapping,
                                fe,
                                support_quadrature,
                                interpolate_update_flags)
        , fe_face_values(mapping, fe, face_quadrature, face_update_flags)
      {}

      // Copy constructor
      ScratchDataSlopeLimiter(const ScratchDataSlopeLimiter<dim> &scratch_data)
        : fe_values(scratch_data.fe_values.get_mapping(),
                    scratch_data.fe_values.get_fe(),
                    scratch_data.fe_values.get_quadrature(),
                    scratch_data.fe_values.get_update_flags())
        , fe_values_neighbor(scratch_data.fe_values_neighbor.get_mapping(),
                             scratch_data.fe_values_neighbor.get_fe(),
                             scratch_data.fe_values_neighbor.get_quadrature(),
                             scratch_data.fe_values_neighbor.get_update_flags())
        , fe_values_interpolate(
            scratch_data.fe_values_interpolate.get_mapping(),
            scratch_data.fe_values_interpolate.get_fe(),
            scratch_data.fe_values_interpolate.get_quadrature(),
            scratch_data.fe_values_interpolate.get_update_flags())
        , fe_face_values(scratch_data.fe_face_values.get_mapping(),
                         scratch_data.fe_face_values.get_fe(),
                         scratch_data.fe_face_values.get_quadrature(),
                         scratch_data.fe_face_values.get_update_flags())
      {}

      FEValues<dim>     fe_values;
      FEValues<dim>     fe_values_neighbor;
      FEValues<dim>     fe_values_interpolate;
      FEFaceValues<dim> fe_face_values;
    };

    struct CopyDataSlopeLimiter
    {
      std::vector<double>                  cell_vector;
      std::vector<types::global_dof_index> local_dof_indices;
      std::vector<types::global_dof_index> local_dof_indices_neighbor;

      template <typename Iterator>
      void
      reinit(const Iterator &cell, unsigned int dofs_per_cell)
      {
        cell_vector.clear();
        cell_vector.resize(dofs_per_cell);

        local_dof_indices.resize(dofs_per_cell);
        local_dof_indices_neighbor.resize(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);
      }
    };

  } // namespace Hydro
} // namespace Sapphire

template <int dim>
Sapphire::Hydro::BurgersEq<dim>::BurgersEq(
  Function<dim>           *initial_condition,
  Function<dim>           *boundary_values,
  Function<dim>           *exact_solution,
  const ParameterParser   &prm,
  const OutputModule<dim> &output_module,
  const double             beta)
  : initial_condition(initial_condition)
  , boundary_values(boundary_values)
  , exact_solution(exact_solution)
  , hd_solver_control(prm)
  , flux(prm, beta)
  , output_module(output_module)
  , beta(beta)
  , mpi_communicator(MPI_COMM_WORLD)
  , mapping()
  , fe(hd_solver_control.fe_degree)
  , dof_handler(triangulation)
  , quadrature_formula(fe.tensor_degree() + 1)
  , face_quadrature_formula(fe.tensor_degree() + 1)
  , pcout(std::cout,
          (Utilities::MPI::this_mpi_process(mpi_communicator) == 0) &&
            (DEBUG_LEVEL >= 0))
  , computing_timer(mpi_communicator,
                    pcout,
                    TimerOutput::never,
                    TimerOutput::wall_times)
{
  AssertDimension(dim, 1);
  DEBUG_PRINT(pcout, 0, "Create BurgersEq");
}

template <int dim>
void
Sapphire::Hydro::BurgersEq<dim>::make_grid()
{
  TimerOutput::Scope t(computing_timer, "Make grid");
  DEBUG_PRINT(pcout, 1, "Make grid");

  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(hd_solver_control.refinement_level);
  DEBUG_PRINT(pcout,
              3,
              "Number of active cells:\t" << triangulation.n_active_cells());
}

template <int dim>
void
Sapphire::Hydro::BurgersEq<dim>::setup_system()
{
  TimerOutput::Scope t(computing_timer, "Setup system");
  DEBUG_PRINT(pcout, 1, "Setup system");

  dof_handler.distribute_dofs(fe);

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
  sparcity_pattern.copy_from(dsp);

  mass_matrix.reinit(sparcity_pattern);
  system_matrix.reinit(sparcity_pattern);

  solution.reinit(dof_handler.n_dofs());
  current_solution.reinit(dof_handler.n_dofs());
  old_solution.reinit(dof_handler.n_dofs());
  dg_vector.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  VectorTools::interpolate(mapping, dof_handler, *initial_condition, solution);

  constraints.clear();
  constraints.close();

  mark_cell_for_limiter.reinit(triangulation.n_active_cells());
  if (hd_solver_control.limiter_criterion == SlopeLimiterCriterion::Always)
    {
      mark_cell_for_limiter = 1.0;
    }
}

template <int dim>
void
Sapphire::Hydro::BurgersEq<dim>::assemble_mass_matrix()
{
  TimerOutput::Scope t(computing_timer, "Assemble mass matrix");
  DEBUG_PRINT(pcout, 1, "Assemble mass matrix");

  MatrixCreator::create_mass_matrix(mapping,
                                    dof_handler,
                                    quadrature_formula,
                                    mass_matrix);
}

template <int dim>
void
Sapphire::Hydro::BurgersEq<dim>::assemble_dg_vector()
{
  TimerOutput::Scope t(computing_timer, "Assemble DG vector");
  DEBUG_PRINT(pcout, 2, "Assemble DG vector");

  dg_vector = 0;
  boundary_values->set_time(current_time);

  using Iterator = typename DoFHandler<dim>::active_cell_iterator;

  const auto cell_worker = [&](const Iterator     &cell,
                               ScratchDataDG<dim> &scratch_data,
                               CopyDataDG         &copy_data) {
    FEValues<dim> &fe_values = scratch_data.fe_values;

    fe_values.reinit(cell);
    const unsigned int n_dofs = fe_values.get_fe().n_dofs_per_cell();

    Tensor<1, dim>      flux_value;
    std::vector<double> current_solution_values(n_dofs);
    fe_values.get_function_values(current_solution, current_solution_values);

    copy_data.reinit(cell, n_dofs);

    for (const unsigned int q_index : fe_values.quadrature_point_indices())
      {
        flux.flux(current_solution_values[q_index], flux_value);

        for (const unsigned int i : fe_values.dof_indices())
          {
            copy_data.cell_vector(i) -= flux_value *
                                        fe_values.shape_grad(i, q_index) *
                                        fe_values.JxW(q_index);
          }
      }
  };

  const auto boundary_worker = [&](const Iterator     &cell,
                                   const unsigned int &face_no,
                                   ScratchDataDG<dim> &scratch_data,
                                   CopyDataDG         &copy_data) {
    scratch_data.fe_face_values.reinit(cell, face_no);
    const FEFaceValuesBase<dim> &fe_face_values = scratch_data.fe_face_values;

    const unsigned int n_dofs = fe_face_values.get_fe().n_dofs_per_cell();

    double              boundary_value;
    std::vector<double> current_solution_values(n_dofs);
    fe_face_values.get_function_values(current_solution,
                                       current_solution_values);

    for (const unsigned int q_index : fe_face_values.quadrature_point_indices())
      {
        const double f_dot_n = beta * current_solution_values[q_index] *
                               current_solution_values[q_index] *
                               fe_face_values.normal_vector(q_index)[0];

        if (f_dot_n > 0.0)
          { // outflow boundary
            for (const unsigned int i : fe_face_values.dof_indices())
              {
                copy_data.cell_vector(i) +=
                  f_dot_n * fe_face_values.shape_value(i, q_index) *
                  fe_face_values.JxW(q_index);
              }
          }
        else
          { // inflow boundary
            for (const unsigned int i : fe_face_values.dof_indices())
              {
                boundary_value = boundary_values->value(
                  fe_face_values.quadrature_point(q_index));
                const double boundary_f_dot_n =
                  beta * boundary_value * boundary_value *
                  fe_face_values.normal_vector(q_index)[0];

                copy_data.cell_vector(i) +=
                  boundary_f_dot_n * fe_face_values.shape_value(i, q_index) *
                  fe_face_values.JxW(q_index);
              }
          }
      }
  };

  const auto face_worker = [&](const Iterator     &cell,
                               const unsigned int &face_no,
                               const unsigned int &subface_no,
                               const Iterator     &neighbor_cell,
                               const unsigned int &neighbor_face_no,
                               const unsigned int &neighbor_subface_no,
                               ScratchDataDG<dim> &scratch_data,
                               CopyDataDG         &copy_data) {
    // supress unused variable warning
    (void)subface_no;
    (void)neighbor_subface_no;

    FEFaceValues<dim> &fe_face_values = scratch_data.fe_face_values;
    fe_face_values.reinit(cell, face_no);

    FEFaceValues<dim> &fe_face_values_neighbor =
      scratch_data.fe_face_values_neighbor;
    fe_face_values_neighbor.reinit(neighbor_cell, neighbor_face_no);

    copy_data.face_data.emplace_back();
    CopyDataFaceDG &copy_data_face = copy_data.face_data.back();

    const unsigned int n_dofs = fe_face_values.get_fe().n_dofs_per_cell();
    copy_data_face.reinit(cell, neighbor_cell, n_dofs);

    std::vector<double> current_solution_values_1(n_dofs);
    fe_face_values.get_function_values(current_solution,
                                       current_solution_values_1);
    std::vector<double> current_solution_values_2(n_dofs);
    fe_face_values_neighbor.get_function_values(current_solution,
                                                current_solution_values_2);

    double         flux_dot_n;

    for (const unsigned int q_index : fe_face_values.quadrature_point_indices())
      {
        flux_dot_n = flux.numerical_flux(current_solution_values_1[q_index],
                                         current_solution_values_2[q_index],
                                         fe_face_values.normal_vector(q_index));

        for (const unsigned int i : fe_face_values.dof_indices())
          {
            copy_data_face.cell_vector_1(i) +=
              flux_dot_n * fe_face_values.shape_value(i, q_index) *
              fe_face_values.JxW(q_index);

            copy_data_face.cell_vector_2(i) -=
              flux_dot_n * fe_face_values_neighbor.shape_value(i, q_index) *
              fe_face_values.JxW(q_index);
          }
      }
  };

  const auto copier = [&](const CopyDataDG &c) {
    constraints.distribute_local_to_global(c.cell_vector,
                                           c.local_dof_indices,
                                           dg_vector);

    for (auto &cdf : c.face_data)
      {
        constraints.distribute_local_to_global(cdf.cell_vector_1,
                                               cdf.local_dof_indices,
                                               dg_vector);
        constraints.distribute_local_to_global(cdf.cell_vector_2,
                                               cdf.local_dof_indices_neighbor,
                                               dg_vector);
      }
  };

  ScratchDataDG<dim> scratch_data(mapping,
                                  fe,
                                  quadrature_formula,
                                  face_quadrature_formula);
  CopyDataDG         copy_data;

  const auto filtered_iterator_range =
    dof_handler.active_cell_iterators() | IteratorFilters::LocallyOwnedCell();
  MeshWorker::mesh_loop(filtered_iterator_range,
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces |
                          MeshWorker::assemble_own_interior_faces_once,
                        boundary_worker,
                        face_worker);
}

template <int dim>
void
Sapphire::Hydro::BurgersEq<dim>::assemble_system()
{
  TimerOutput::Scope t(computing_timer, "Assemble system");
  DEBUG_PRINT(pcout, 2, "Assemble system");

  system_matrix = 0;
  system_rhs    = 0;

  /** Nothing to do here, RHS of equation is zero */
}

template <int dim>
void
Sapphire::Hydro::BurgersEq<dim>::slope_limiter()
{
  if (hd_solver_control.limiter_criterion == SlopeLimiterCriterion::Never)
    return;
  TimerOutput::Scope t(computing_timer, "Slope limiter");
  DEBUG_PRINT(pcout, 2, "Slope limiter");

  Vector<double> limited_solution(dof_handler.n_dofs());

  // boundary_values->set_time(current_time);

  using Iterator = typename DoFHandler<dim>::active_cell_iterator;

  const auto cell_worker = [&](const Iterator               &cell,
                               ScratchDataSlopeLimiter<dim> &scratch_data,
                               CopyDataSlopeLimiter         &copy_data) {
    FEValues<dim> &fe_values             = scratch_data.fe_values;
    FEValues<dim> &fe_values_neighbor    = scratch_data.fe_values_neighbor;
    FEValues<dim> &fe_values_interpolate = scratch_data.fe_values_interpolate;
    FEFaceValues<dim> &fe_face_values    = scratch_data.fe_face_values;

    fe_values.reinit(cell);
    fe_values_interpolate.reinit(cell);

    const unsigned int n_dofs = fe_values.get_fe().n_dofs_per_cell();
    copy_data.reinit(cell, n_dofs);

    std::vector<types::global_dof_index> &local_dof_indices =
      copy_data.local_dof_indices;
    std::vector<types::global_dof_index> &local_dof_indices_neighbor =
      copy_data.local_dof_indices_neighbor;
    const Point<dim> cell_center = cell->center();

    bool limit_cell = true;
    if (hd_solver_control.limiter_criterion == SlopeLimiterCriterion::Never)
      {
        limit_cell = false;
      }
    if (!limit_cell)
      {
        mark_cell_for_limiter[cell->active_cell_index()] = 0.0;
        for (const unsigned int i : fe_values.dof_indices())
          {
            copy_data.cell_vector[i] = solution[local_dof_indices[i]];
          }
        return;
      }

    // Calculate the average of current cell
    double         cell_average = 0;
    Tensor<1, dim> cell_average_grad;
    for (const unsigned int q_index : fe_values.quadrature_point_indices())
      {
        for (const unsigned int i : fe_values.dof_indices())
          {
            cell_average += solution[local_dof_indices[i]] *
                            fe_values.shape_value(i, q_index) *
                            fe_values.JxW(q_index);
            cell_average_grad += solution[local_dof_indices[i]] *
                                 fe_values.shape_grad(i, q_index) *
                                 fe_values.JxW(q_index);
          }
      }
    cell_average /= cell->measure();
    cell_average_grad /= cell->measure();

    // Calculate the averages of the neighbour cells
    std::vector<double>         neighbor_cell_averages(cell->n_faces());
    std::vector<Tensor<1, dim>> neighbor_distance(cell->n_faces());
    unsigned int                n_neighbors = 0;
    for (const auto face_no : cell->face_indices())
      {
        if (!cell->at_boundary(face_no) && cell->neighbor(face_no)->is_active())
          {
            auto neighbor = cell->neighbor(face_no);
            neighbor->get_dof_indices(local_dof_indices_neighbor);
            fe_values_neighbor.reinit(neighbor);

            neighbor_distance[n_neighbors] =
              neighbor->center() - cell->center();

            neighbor_cell_averages[n_neighbors] = 0;
            for (const unsigned int q_index :
                 fe_values_neighbor.quadrature_point_indices())
              {
                for (const unsigned int i : fe_values_neighbor.dof_indices())
                  {
                    neighbor_cell_averages[n_neighbors] +=
                      solution[local_dof_indices_neighbor[i]] *
                      fe_values_neighbor.shape_value(i, q_index) *
                      fe_values_neighbor.JxW(q_index);
                  }
              }
            neighbor_cell_averages[n_neighbors] /= neighbor->measure();
            n_neighbors++;
          }
      }

    // Check if the cell is meets the criteria for limiting
    if (hd_solver_control.limiter_criterion ==
        SlopeLimiterCriterion::GerneralizedSlopeLimiter)
      {
        limit_cell = false;

        double              face_flux  = 0;
        double              face_value = 0;
        std::vector<double> tmp_minmod(n_neighbors + 1);

        unsigned int i_neighbor = 0;
        for (const auto face_no : cell->face_indices())
          {
            if (!cell->at_boundary(face_no) &&
                cell->neighbor(face_no)->is_active())
              {
                fe_face_values.reinit(cell, face_no);

                Tensor<1, dim> direction = neighbor_distance[i_neighbor] /
                                           neighbor_distance[i_neighbor].norm();

                // Calculate the average of current face
                double weight = 0;
                face_value    = 0;
                for (const unsigned int q_index :
                     fe_face_values.quadrature_point_indices())
                  {
                    for (const unsigned int i : fe_face_values.dof_indices())
                      {
                        face_value += solution[local_dof_indices[i]] *
                                      fe_face_values.shape_value(i, q_index) *
                                      fe_face_values.JxW(q_index);
                      }
                    weight += fe_face_values.JxW(q_index);
                  }
                face_value /= weight;

                tmp_minmod[0] = (face_value - cell_average);
                for (unsigned int i2 = 0; i2 < n_neighbors; i2++)
                  {
                    tmp_minmod[i2 + 1] =
                      (neighbor_cell_averages[i2] - cell_average) *
                      (neighbor_distance[i2] * direction) /
                      neighbor_distance[i2].norm();
                  }
                face_flux = cell_average + minmod(tmp_minmod);

                if (std::abs(face_value - face_flux) > 1e-10)
                  {
                    limit_cell = true;
                    break;
                  }
                i_neighbor++;
              }
          }
        if (!limit_cell)
          AssertDimension(i_neighbor, n_neighbors);
      }

    mark_cell_for_limiter[cell->active_cell_index()] = limit_cell ? 1.0 : 0.0;

    if (limit_cell)
      {
        // Calculate the limited slope
        Tensor<1, dim> limited_slope;
        compute_limited_slope(cell_average,
                              cell_average_grad,
                              neighbor_cell_averages,
                              neighbor_distance,
                              n_neighbors,
                              limited_slope,
                              hd_solver_control);

        // To calculate the updated dof-values, we use a similar functionality
        // as VectroTools::interpolate

        // Calculate the support point values
        const std::vector<Point<dim>> &support_points =
          fe_values_interpolate.get_present_fe_values().get_quadrature_points();
        AssertDimension(support_points.size(), n_dofs);
        std::vector<Vector<double>> support_point_values(n_dofs,
                                                         Vector<double>(1));

        for (unsigned int i = 0; i < n_dofs; ++i)
          {
            support_point_values[i] =
              cell_average + (support_points[i] - cell_center) * limited_slope;
          }

        // Transformation in these cases needed, not implemented!
        Assert(fe.conforming_space != FiniteElementData<dim>::Hcurl &&
                 fe.conforming_space != FiniteElementData<dim>::Hdiv &&
                 fe.conforming_space != FiniteElementData<dim>::H1,
               ExcNotImplemented());

        fe.convert_generalized_support_point_values_to_dof_values(
          support_point_values, copy_data.cell_vector);
      }
    else
      {
        // No limiting, just copy the values
        for (const unsigned int i : fe_values.dof_indices())
          {
            copy_data.cell_vector[i] = solution[local_dof_indices[i]];
          }
      }
  };

  const auto copier = [&](const CopyDataSlopeLimiter &c) {
    constraints.distribute_local_to_global(c.cell_vector,
                                           c.local_dof_indices,
                                           limited_solution);
  };

  Quadrature<dim> support_quadrature(fe.get_generalized_support_points());

  ScratchDataSlopeLimiter<dim> scratch_data(mapping,
                                            fe,
                                            quadrature_formula,
                                            support_quadrature,
                                            face_quadrature_formula);
  CopyDataSlopeLimiter         copy_data;

  const auto filtered_iterator_range =
    dof_handler.active_cell_iterators() | IteratorFilters::LocallyOwnedCell();
  MeshWorker::mesh_loop(filtered_iterator_range,
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells);

  solution = limited_solution;
}

template <int dim>
void
Sapphire::Hydro::BurgersEq<dim>::perform_time_step()
{
  TimerOutput::Scope t(computing_timer, "Time step");
  DEBUG_PRINT(pcout, 2, "Time step");

  old_solution     = solution;
  current_solution = old_solution;
  Vector<double> tmp(dof_handler.n_dofs());

  double time_step = hd_solver_control.time_step;

  switch (hd_solver_control.scheme)
    {
      case TimeSteppingScheme::ForwardEuler:
        {
          assemble_system();
          assemble_dg_vector();

          system_rhs.add(-1.0, dg_vector);
          system_rhs *= time_step;

          mass_matrix.vmult(tmp, current_solution);
          system_rhs.add(1.0, tmp);

          system_matrix.copy_from(mass_matrix);

          solve_linear_system();

          slope_limiter();

          break;
        }

      case TimeSteppingScheme::ExplicitRK:
        {
          //  Butcher's array
          Vector<double> a({0.5, 0.5, 1.});
          Vector<double> b({1. / 6, 1. / 3, 1. / 3, 1. / 6});
          Vector<double> c({0., 0.5, 0.5, 1.});
          int            i = 0;

          Vector<double> k1(dof_handler.n_dofs());
          i                = 0;
          tmp              = 0;
          current_time     = time + c[i] * time_step;
          current_solution = old_solution;
          assemble_system();
          assemble_dg_vector();
          system_matrix.copy_from(mass_matrix);
          system_rhs.add(-1.0, dg_vector);
          solve_linear_system();
          slope_limiter();
          k1 = solution;

          Vector<double> k2(dof_handler.n_dofs());
          i                = 1;
          tmp              = 0;
          current_time     = time + c[i] * time_step;
          current_solution = old_solution;
          current_solution.add(a[i - 1] * time_step, k1);
          assemble_system();
          assemble_dg_vector();
          system_matrix.copy_from(mass_matrix);
          system_rhs.add(-1.0, dg_vector);
          solve_linear_system();
          slope_limiter();
          k2 = solution;

          Vector<double> k3(dof_handler.n_dofs());
          i                = 2;
          tmp              = 0;
          current_time     = time + c[i] * time_step;
          current_solution = old_solution;
          current_solution.add(a[i - 1] * time_step, k2);
          assemble_system();
          assemble_dg_vector();
          system_matrix.copy_from(mass_matrix);
          system_rhs.add(-1.0, dg_vector);
          solve_linear_system();
          slope_limiter();
          k3 = solution;

          Vector<double> k4(dof_handler.n_dofs());
          i                = 3;
          tmp              = 0;
          current_time     = time + c[i] * time_step;
          current_solution = old_solution;
          current_solution.add(a[i - 1] * time_step, k3);
          assemble_system();
          assemble_dg_vector();
          system_matrix.copy_from(mass_matrix);
          system_rhs.add(-1.0, dg_vector);
          solve_linear_system();
          slope_limiter();
          k4 = solution;

          solution = old_solution;
          solution.add(b[0] * time_step, k1);
          solution.add(b[1] * time_step, k2);
          solution.add(b[2] * time_step, k3);
          solution.add(b[3] * time_step, k4);
          slope_limiter();
          break;
        }

      default:
        Assert(false, ExcNotImplemented());
        break;
    }

  time += time_step;
}

template <int dim>
void
Sapphire::Hydro::BurgersEq<dim>::solve_linear_system()
{
  TimerOutput::Scope t(computing_timer, "Solve linear system");
  DEBUG_PRINT(pcout, 2, "Solve linear system");

  SolverControl solver_control(hd_solver_control.max_iterations,
                               hd_solver_control.tolerance);
  // SolverCG<Vector<double>> solver(solver_control);
  SolverRichardson<Vector<double>> solver(solver_control);

  // PreconditionBlockJacobi preconditioner;
  // preconditioner.initialize(system_matrix);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix);

  solver.solve(system_matrix, solution, system_rhs, preconditioner);

  // constraints.distribute(solution);

  DEBUG_PRINT(pcout,
              3,
              "Solver converged in " << solver_control.last_step()
                                     << " iterations.");
}

template <int dim>
void
Sapphire::Hydro::BurgersEq<dim>::output_results()
{
  DEBUG_PRINT(pcout, 2, "Output results");

  Vector<double> exact_solution_values(dof_handler.n_dofs());
  exact_solution->set_time(time);
  VectorTools::interpolate(mapping,
                           dof_handler,
                           *exact_solution,
                           exact_solution_values);

  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "Solution");
  data_out.add_data_vector(exact_solution_values, "ExactSolution");
  data_out.add_data_vector(mark_cell_for_limiter, "LimitCell");

  data_out.build_patches(fe.tensor_degree());

  output_module.write_results(data_out, timestep_number);
}

template <int dim>
void
Sapphire::Hydro::BurgersEq<dim>::process_results()
{
  TimerOutput::Scope t(computing_timer, "Process results");
  DEBUG_PRINT(pcout, 2, "Process results");

  Vector<float> difference_per_cell(triangulation.n_active_cells());
  exact_solution->set_time(time);

  // Use different quadrature for error computation
  const QTrapezoid<1>  q_trapez;
  const QIterated<dim> q_iterated(q_trapez, fe.degree * 2 + 1);

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    *exact_solution,
                                    difference_per_cell,
                                    q_iterated,
                                    VectorTools::L2_norm);
  float L2_error = VectorTools::compute_global_error(triangulation,
                                                     difference_per_cell,
                                                     VectorTools::L2_norm);
  DEBUG_PRINT(pcout, 3, "L2 error:\t\t" << L2_error);

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    *exact_solution,
                                    difference_per_cell,
                                    q_iterated,
                                    VectorTools::Linfty_norm);
  float Linf_error =
    VectorTools::compute_global_error(triangulation,
                                      difference_per_cell,
                                      VectorTools::Linfty_norm);
  DEBUG_PRINT(pcout, 3, "L-infinity error:\t" << Linf_error);

  error_with_time.push_back(L2_error);
}

template <int dim>
void
Sapphire::Hydro::BurgersEq<dim>::init()
{
  DEBUG_PRINT(pcout, 0, "Init BurgersEq");
  time            = 0.0;
  timestep_number = 0;
  error_with_time.clear();
  error_with_time.reserve(
    (unsigned int)(hd_solver_control.end_time / hd_solver_control.time_step) /
    output_module.output_frequency);

  make_grid();
  setup_system();
  assemble_mass_matrix();
  // assemble_system();
  {
    TimerOutput::Scope t(computing_timer, "Output results");
    output_results();
  }
}

template <int dim>
void
Sapphire::Hydro::BurgersEq<dim>::do_timestep()
{
  DEBUG_PRINT(pcout,
              1,
              "Timestep " << timestep_number + 1
                          << " (time = " << time + hd_solver_control.time_step
                          << "/" << hd_solver_control.end_time << ")");
  perform_time_step();
  timestep_number++;
  if (timestep_number % output_module.output_frequency == 0)
    {
      {
        TimerOutput::Scope t(computing_timer, "Output results");
        output_results();
      }
      process_results();
    }
}

template <int dim>
void
Sapphire::Hydro::BurgersEq<dim>::run()
{
  DEBUG_PRINT(pcout, 0, "Run BurgersEq");
  init();
  while (time < hd_solver_control.end_time)
    {
      do_timestep();
    }

  computing_timer.print_summary();
  computing_timer.reset();

  if (pcout.is_active() && DEBUG_LEVEL >= 2)
    {
      pcout << "   L2 error with time:" << std::endl;
      pcout << "   ";
      for (const auto &e : error_with_time)
        {
          // pcout << e << " ";
          printf("%.2e ", e);
        }
      pcout << std::endl;
    }
}

// explicit instantiation
template class Sapphire::Hydro::BurgersEq<1>;
