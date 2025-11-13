"""Plot steady-state-parallel-shock-scaled example."""

import paraview.simple as ps
from sapphireppplot import vfp, pvplot, transform


def main() -> dict:
    """Plot steady-state-parallel-shock-scaled example."""
    plot_properties = vfp.PlotPropertiesVFP(
        dimension=2,
        momentum=True,
        scaled_distribution_function=True,
        lms_indices=[[0, 0, 0]],
        axes_stretch=[1.0, 2e3, 1.0],
        legend_location="BottomLeft",
    )

    results_folder, prm, solution, animation_scene = vfp.load_solution(
        plot_properties,
        path_prefix="$SAPPHIREPP_RESULTS/steady-state-parallel-shock",
    )

    # region Physical parameters for the analytical solution
    source = float(prm["Physical parameters"]["Q"])  # type: ignore
    p_inj = float(prm["Physical parameters"]["p_inj"])  # type: ignore
    r = float(prm["Physical parameters"]["compression ratio"])  # type: ignore
    u_one = float(prm["Physical parameters"]["u_sh"])  # type: ignore
    bounds = solution.GetDataInformation().GetBounds()
    # endregion

    # region Plot 2D render view
    layout_2d, render_view_2d = vfp.plot_f_lms_2d(
        solution,
        results_folder,
        "steady-state-parallel-shock-scaled-2D",
        plot_properties,
        value_range=[1e-6, 60.0],
    )
    # endregion

    # region Scale distribution function to p^4 f
    solution_scaled, plot_properties_scaled = vfp.scale_distribution_function(
        solution,
        plot_properties.replace(interpol=True),
        lms_indices=[[0, 0, 0]],
        spectral_index=4,
    )
    # endregion

    # region Plot f(p)
    scaled_f_000_ana_p = ps.PythonCalculator(
        registrationName="scaled_f_000_ana_p", Input=solution_scaled
    )
    scaled_f_000_ana_p.ArrayName = "interpol_p^s f_000"
    scaled_f_000_ana_p.UseMultilineExpression = 1
    scaled_f_000_ana_p.MultilineExpression = f"""
    # Physical Parameters
    source = {source}
    p_inj = {p_inj}
    r = {r}
    u_one = {u_one}

    N = 3*source/(numpy.sqrt(4 * numpy.pi) * u_one * p_inj**3) \
        * r/(r - 1) * (p_inj)**(3*r/(r-1))

    outputArray = N * numpy.ones(inputs[0].PointData['g_000'].shape[0])

    return outputArray
    """

    plot_over_line_p = transform.plot_over_line(
        scaled_f_000_ana_p,
        direction="y",
        x_range=[0.0, bounds[3]],
        offset=[0.1, 0.0, 0.0],
        results_folder=results_folder,
        filename="scaled-particle-spectrum",
        plot_properties=plot_properties_scaled,
    )
    layout_p = ps.CreateLayout("scaled-particle-spectrum")
    line_chart_view_p = pvplot.plot_line_chart_view(
        plot_over_line_p,
        layout_p,
        x_label=r"$\ln p$",
        y_label=r"$p^s f_{000}$",
        x_array_name="Points_Y",
        visible_lines=["g_000", "p^s f_000", "interpol_p^s f_000"],
        log_y_scale=False,
        plot_properties=plot_properties_scaled,
    )
    pvplot.save_screenshot(
        layout_p,
        results_folder,
        "scaled-particle-spectrum",
        plot_properties_scaled,
    )
    # endregion

    return locals()


if __name__ in ["__main__", "__vtkconsole__"]:
    results = main()
    # Make all results available as global variables in a vtkconsole
    globals().update(results)
