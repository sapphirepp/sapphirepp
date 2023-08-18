#include "postprocessor.h"

#include "flux.h"

template <int dim>
Sapphire::Hydro::Postprocessor<dim>::Postprocessor()
  : do_schlieren_plot(dim == 2)
{}



template <int dim>
void
Sapphire::Hydro::Postprocessor<dim>::evaluate_vector_field(
  const dealii::DataPostprocessorInputs::Vector<dim> &inputs,
  std::vector<Vector<double>>                        &computed_quantities) const
{
  const unsigned int n_evaluation_points = inputs.solution_values.size();

  if (do_schlieren_plot == true)
    Assert(inputs.solution_gradients.size() == n_evaluation_points,
           ExcInternalError());

  Assert(computed_quantities.size() == n_evaluation_points, ExcInternalError());
  Assert(inputs.solution_values[0].size() == dim + 2, ExcInternalError());
  Assert(computed_quantities[0].size() ==
           dim + 2 + (do_schlieren_plot == true ? 1 : 0),
         ExcInternalError());

  for (unsigned int p = 0; p < n_evaluation_points; ++p)
    {
      Tensor<1, dim + 2> solution;
      for (unsigned int d = 0; d < dim + 2; ++d)
        solution[d] = inputs.solution_values[p](d);

      const double         density  = solution[0];
      const Tensor<1, dim> velocity = euler_velocity<dim>(solution);
      const double         pressure = euler_pressure<dim>(solution);

      for (unsigned int d = 0; d < dim; ++d)
        computed_quantities[p](d) = velocity[d];
      computed_quantities[p](dim)     = pressure;
      computed_quantities[p](dim + 1) = std::sqrt(gamma * pressure / density);

      if (do_schlieren_plot == true)
        computed_quantities[p](dim + 2) =
          inputs.solution_gradients[p][0] * inputs.solution_gradients[p][0];
    }
}



template <int dim>
std::vector<std::string>
Sapphire::Hydro::Postprocessor<dim>::get_names() const
{
  std::vector<std::string> names;
  for (unsigned int d = 0; d < dim; ++d)
    names.emplace_back("velocity");
  names.emplace_back("pressure");
  names.emplace_back("speed_of_sound");

  if (do_schlieren_plot == true)
    names.emplace_back("schlieren_plot");

  return names;
}



template <int dim>
std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
Sapphire::Hydro::Postprocessor<dim>::get_data_component_interpretation() const
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation;
  for (unsigned int d = 0; d < dim; ++d)
    interpretation.push_back(
      DataComponentInterpretation::component_is_part_of_vector);
  interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  interpretation.push_back(DataComponentInterpretation::component_is_scalar);

  if (do_schlieren_plot == true)
    interpretation.push_back(DataComponentInterpretation::component_is_scalar);

  return interpretation;
}



template <int dim>
dealii::UpdateFlags
Sapphire::Hydro::Postprocessor<dim>::get_needed_update_flags() const
{
  if (do_schlieren_plot == true)
    return update_values | update_gradients;
  else
    return update_values;
}

// explicit instantiation
template class Sapphire::Hydro::Postprocessor<1>;
template class Sapphire::Hydro::Postprocessor<2>;
template class Sapphire::Hydro::Postprocessor<3>;