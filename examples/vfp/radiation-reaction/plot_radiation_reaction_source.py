"""Plot radiation-reaction example."""

from sapphireppplot import vfp, utils


def main() -> dict:
    """Plot radiation-reaction example."""
    plot_properties = vfp.PlotPropertiesVFP(
        dimension=1,
        momentum=True,
        lms_indices=[[0, 0, 0]],
        scaled_distribution_function=True,
        line_colors={
            "g_000": utils.sapphirepp_colors()[0],
            "numeric_g_000": utils.sapphirepp_colors()[0],
            "interpol_g_000": utils.sapphirepp_colors()[0],
        },
        # For test runs with analytic solution:
        prefix_numeric=True,
        interpol=True,
    )
    # plot_properties.convert_lnp_to_p()

    results_folder, prm, solution, animation_scene = vfp.load_solution(
        plot_properties,
        results_folder="$SAPPHIREPP_RESULTS/radiation-reaction-with-source/",
    )

    solution_scaled, plot_properties_scaled = vfp.scale_distribution_function(
        solution, plot_properties
    )

    layout, line_chart_view = vfp.plot_f_lms_1d(
        solution_scaled,
        results_folder,
        "radiation-reaction-with-source",
        plot_properties_scaled,
        value_range=[0, 0.5],
        # log_y_scale=True,
        save_animation=True,
    )

    return locals()


if __name__ in ["__main__", "__vtkconsole__"]:
    results = main()
    # Make all results available as global variables in a vtkconsole
    globals().update(results)
