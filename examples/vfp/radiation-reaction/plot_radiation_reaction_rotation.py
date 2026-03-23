"""Plot radiation-reaction example with rotation."""

from sapphireppplot import vfp, utils


def main() -> dict:
    """Plot radiation-reaction example with rotation."""
    plot_properties = vfp.PlotPropertiesVFP(
        dimension=1,
        momentum=True,
        lms_indices=[[1, 1, 0], [1, 1, 1]],
        scaled_distribution_function=True,
        line_colors={
            "g_110": utils.sapphirepp_colors()[0],
            "g_111": utils.sapphirepp_colors()[1],
            "numeric_g_110": utils.sapphirepp_colors()[0],
            "numeric_g_111": utils.sapphirepp_colors()[1],
            "project_g_110": utils.sapphirepp_colors()[0],
            "project_g_111": utils.sapphirepp_colors()[1],
        },
        # For test runs with analytic solution:
        prefix_numeric=True,
        project=True,
    )
    plot_properties.convert_lnp_to_p()

    results_folder, prm, solution, animation_scene = vfp.load_solution(
        plot_properties,
        results_folder="$SAPPHIREPP_RESULTS/radiation-reaction-with-rotation/",
    )

    solution_scaled, plot_properties_scaled = vfp.scale_distribution_function(
        solution, plot_properties
    )

    layout, line_chart_view = vfp.plot_f_lms_1d(
        solution_scaled,
        results_folder,
        "radiation-reaction-with-rotation",
        plot_properties_scaled,
        value_range=[-3.2, 3.2],
        save_animation=True,
    )

    return locals()


if __name__ in ["__main__", "__vtkconsole__"]:
    results = main()
    # Make all results available as global variables in a vtkconsole
    globals().update(results)
