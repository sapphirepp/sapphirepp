/**
 * @file postprocessor.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement Postprocessor class
 * @version 0.1
 * @date 2023-08-18
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef HYDROSOLVER_POSTPROCESSOR_H
#define HYDROSOLVER_POSTPROCESSOR_H

#include <deal.II/numerics/data_out.h>

namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;

    template <int dim>
    class Postprocessor : public DataPostprocessor<dim>
    {
    public:
      Postprocessor();

      virtual void
      evaluate_vector_field(
        const DataPostprocessorInputs::Vector<dim> &inputs,
        std::vector<Vector<double>> &computed_quantities) const override;

      virtual std::vector<std::string>
      get_names() const override;

      virtual std::vector<
        DataComponentInterpretation::DataComponentInterpretation>
      get_data_component_interpretation() const override;

      virtual UpdateFlags
      get_needed_update_flags() const override;

    private:
      const bool do_schlieren_plot;
    };
  } // namespace Hydro
} // namespace Sapphire

#endif