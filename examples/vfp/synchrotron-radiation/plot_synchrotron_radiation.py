"""Plot synchrotron-radiation example."""

import numpy as np
import paraview.simple as ps
from sapphireppplot import vfp, pvplot, transform


def main() -> dict:
    """Plot synchrotron-radiation example."""
    plot_properties = vfp.PlotPropertiesVFP(
        dimension=2,
        momentum=True,
        scaled_distribution_function=True,
        lms_indices=[[0, 0, 0]],
        axes_stretch=[1.0, 2e3, 1.0],
        # legend_location="BottomLeft",
        sampling_resolution=1e-8,
    )

    results_folder, prm, solution, animation_scene = vfp.load_solution(
        plot_properties,
        path_prefix="$SAPPHIREPP_RESULTS/synchrotron-radiation",
        results_folder="/home/schulze/Documents/PhD/Code/sapphirepp/results/synchrotron-radiation",
    )

    # region Physical parameters for the analytical solution
    source = float(prm["Physical parameters"]["Q"])  # type: ignore
    p_inj = float(prm["Physical parameters"]["p_inj"])  # type: ignore
    r = float(prm["Physical parameters"]["compression ratio"])  # type: ignore
    u_one = float(prm["Physical parameters"]["u_sh"])  # type: ignore
    p_hat = 4 * p_inj  # spatial dependence of f_000 is shown at p = p_hat
    bounds = solution.GetDataInformation().GetBounds()
    # endregion

    # region Plot 2D render view
    layout_2d, render_view_2d = vfp.plot_f_lms_2d(
        solution,
        results_folder,
        "synchrotron-radiation-2D",
        plot_properties,
        value_range=[1e-6, 60.0],
    )
    # endregion

    # region Scale distribution function to p^4 f
    solution_scaled, plot_properties_scaled = vfp.scale_distribution_function(
        solution,
        plot_properties,
        lms_indices=[[0, 0, 0]],
        spectral_index=4,
    )
    # endregion

    # region Plot f(x)
    plot_over_line_x, layout_x, line_chart_view_x = vfp.plot_f_lms_over_x(
        solution_scaled,
        results_folder,
        "spatial-distribution",
        plot_properties_scaled,
        direction="x",
        offset=[0.0, np.log(p_hat), 0.0],
        # x_range=[bounds[0] * 0.1, bounds[1]],  # 10% of the upstream
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

    plot_properties_p = plot_properties_scaled.replace(
        interpol=1,
        annotation_project_interpol="exact",
    )
    plot_over_line_p, layout_p, line_chart_view_p = vfp.plot_f_lms_over_p(
        scaled_f_000_ana_p,
        results_folder,
        "synchrotron-radiation-spectrum",
        plot_properties_p,
        offset=[0.1, 0.0, 0.0],
        value_range=[5e6, 2e7],
        log_y_scale=True,
    )
    # endregion

    return locals()


if __name__ in ["__main__", "__vtkconsole__"]:
    results = main()
    # Make all results available as global variables in a vtkconsole
    globals().update(results)
