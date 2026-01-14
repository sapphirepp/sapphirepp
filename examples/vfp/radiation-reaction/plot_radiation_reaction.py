"""Plot radiation-reaction example."""

import numpy as np
import paraview.simple as ps
from sapphireppplot import vfp, pvplot, transform


def main() -> dict:
    """Plot radiation-reaction example."""
    plot_properties = vfp.PlotPropertiesVFP(
        dimension=1,
        momentum=True,
        lms_indices=[[0, 0, 0], [1, 0, 0]],
        scaled_distribution_function=True,
        prefix_numeric=True,
        project=True,
    )

    results_folder, prm, solution, animation_scene = vfp.load_solution(
        plot_properties,
        # path_prefix="$SAPPHIREPP_RESULTS/parallel-shock",
    )

    layout, line_chart_view = vfp.plot_f_lms_1d(
        solution,
        results_folder,
        "radiation-reaction",
        plot_properties,
        value_range=[1e-6, 1],
        log_x_scale=False,
        log_y_scale=True,
        # save_animation=True,
    )

    # region Plot over time
    ln_p_hat = 16.0
    probe_location, plot_properties_t = transform.probe_location(
        solution, [ln_p_hat, 0.0, 0.0], plot_properties
    )

    plot_over_time, plot_properties_t = transform.plot_over_time(
        probe_location,
        results_folder=results_folder,
        filename="temporal-evolution",
        plot_properties_in=plot_properties_t,
    )

    layout_t = ps.CreateLayout("temporal-evolution")
    line_chart_view_t = pvplot.plot_line_chart_view(
        plot_over_time,
        layout_t,
        x_label=r"$t$",
        y_label=r"$f_{lms}$",
        x_array_name="Time",
        visible_lines=[
            f"{name} (id=0)" for name in plot_properties.series_names
        ],
        log_y_scale=True,
        plot_properties=plot_properties_t,
    )
    pvplot.save_screenshot(
        layout_t, results_folder, "temporal-evolution", plot_properties_t
    )
    # endregion

    return locals()


if __name__ in ["__main__", "__vtkconsole__"]:
    results = main()
    # Make all results available as global variables in a vtkconsole
    globals().update(results)
