"""Plot steady-state-radiation-reaction-shock example."""

import numpy as np
from sapphireppplot import vfp


def main() -> dict:
    """Plot steady-state-radiation-reaction-shock example with p^3 rescaling."""
    plot_properties = vfp.PlotPropertiesVFP(
        dimension=2,
        momentum=True,
        scaled_distribution_function=True,
        axes_stretch=(5e-10, 1, 1),
        sampling_pattern="uniform",
        sampling_resolution=10000,
        # sampling_pattern='center', # avoid oscillation artifacts
        # sampling_resolution=1.e-12, # choose on run-by-run basis
    )

    results_folder, prm, solution, animation_scene = vfp.load_solution(
        plot_properties,
        path_prefix="$SAPPHIREPP_RESULTS/steady-state-radiation-reaction-shock",
    )

    # region Plot 2D render view
    layout_2d, render_view_2d = vfp.plot_f_lms_2d(
        solution,
        results_folder,
        "steady-state-radiation-reaction-shock-2D",
        plot_properties,
        lms_index=(0, 0, 0),
        value_range=(1e-4, 1e1),
    )
    # endregion

    # region Plot spatial distribution g(x)
    plot_over_line_x, layout_x, line_chart_view_x = vfp.plot_f_lms_over_x(
        solution,
        results_folder,
        "steady-state-radiation-reaction-shock-x",
        plot_properties,
        lms_indices=[(0, 0, 0), (1, 0, 0)],
        direction="x",
        offset=(0, 17.0, 0),  # choose momentum cutout
        x_label=r"$x$",
        x_range=(-1e3, 0.5e9),  # choose plotted range
    )
    # endregion

    # region Plot F(p)
    # injection momentum (lower p boundary) as normalization:
    p_inj = np.exp(float(prm["VFP"]["Mesh"]["Point 1"].split(",")[-1]))

    # momentum plot scaled and normalized, p^4 f(p)
    solution_scaled, plot_properties_scaled = vfp.scale_distribution_function(
        solution,
        plot_properties,
        spectral_index=4.0,
        # normalization = 1/p_inj
    )

    plot_over_line_p, layout_p, line_chart_view_p = vfp.plot_f_lms_over_p(
        solution_scaled,
        results_folder,
        "steady-state-radiation-reaction-shock-p",
        plot_properties_scaled,
        lms_indices=[(0, 0, 0)],
        offset=(0.1, 0, 0),
        value_range=(1e6, 8e6),
    )
    # endregion

    return locals()


if __name__ in ["__main__", "__vtkconsole__"]:
    results = main()
    globals().update(results)
