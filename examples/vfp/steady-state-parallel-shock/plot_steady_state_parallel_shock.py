"""Plot steady-state-parallel-shock example."""

import numpy as np
import paraview.simple as ps
from sapphireppplot import vfp


def main() -> dict:
    """Plot steady-state-parallel-shock example."""
    plot_properties = vfp.PlotPropertiesVFP(
        dimension=2,
        momentum=True,
        lms_indices=[[0, 0, 0]],
        axes_stretch=[1.0, 2e3, 1.0],
        sampling_resolution=1e-8,
    )

    results_folder, prm, solution, animation_scene = vfp.load_solution(
        plot_properties,
        path_prefix="$SAPPHIREPP_RESULTS/steady-state-parallel-shock",
    )

    # region Plot 2D render view
    layout_2d, render_view_2d = vfp.plot_f_lms_2d(
        solution,
        results_folder,
        "steady-state-parallel-shock-2D",
        plot_properties,
        value_range=[1e-6, 30.0],
    )
    # endregion

    # region Physical parameters for the analytical solution
    source = float(prm["Physical parameters"]["Q"])  # type: ignore
    p_inj = float(prm["Physical parameters"]["p_inj"])  # type: ignore
    r = float(prm["Physical parameters"]["compression ratio"])  # type: ignore
    u_one = float(prm["Physical parameters"]["u_sh"])  # type: ignore
    nu_0 = float(prm["Physical parameters"]["nu0"])  # type: ignore
    p_hat = 50  # spatial dependence of f_000 is shown at p = p_hat
    bounds = solution.GetDataInformation().GetBounds()
    # endregion

    # region Plot f(p)
    f_000_ana_p = ps.PythonCalculator(
        registrationName="f_000_ana_p", Input=solution
    )
    f_000_ana_p.ArrayName = "interpol_f_000"
    f_000_ana_p.UseMultilineExpression = 1
    f_000_ana_p.MultilineExpression = f"""
    # Physical Parameters
    source = {source}
    p_inj = {p_inj}
    r = {r}
    u_one = {u_one}

    p_coords = numpy.exp(inputs[0].Points[:,1])
    outputArray = 3*source/(numpy.sqrt(4 * numpy.pi) * u_one * p_inj**3) \
        * r/(r - 1) * (p_coords/p_inj)**(-3*r/(r-1))

    return outputArray
    """

    plot_properties_p = plot_properties.replace(interpol=True)
    plot_over_line_p, layout_p, line_chart_view_p = vfp.plot_f_lms_over_p(
        f_000_ana_p,
        results_folder,
        "particle-spectrum",
        plot_properties_p,
        offset=[0.1, 0, 0],
        x_range=[0.0, bounds[3]],
        value_range=[1e-6, 14.5],
    )
    # endregion

    # region Plot f(x)
    f_000_ana_x = ps.PythonCalculator(
        registrationName="f_000_ana_x", Input=solution
    )
    f_000_ana_x.ArrayName = "interpol_f_000"
    f_000_ana_x.UseMultilineExpression = 1
    f_000_ana_x.MultilineExpression = f"""
    # Physical Parameters
    source = {source}
    p_inj = {p_inj}
    r = {r}
    u_one = {u_one}
    nu_0 = {nu_0}
    p_hat = {p_hat}

    # Derived quantities
    v = p_hat/numpy.sqrt(p_hat**2 + 1)   # magnitude of the particle velocity
    nu = nu_0/p_hat                      # scattering frequency for p_hat
    # Normalization of the distribution function
    N = 3*source/(numpy.sqrt(4 * numpy.pi) * u_one * p_inj**3) \
        * r/(r - 1) * (p_hat/p_inj)**(-3*r/(r-1))

    x_coords = inputs[0].Points[:,0]
    outputArray = numpy.zeros(x_coords.shape[0])

    # Upstream
    outputArray[x_coords < 0] = N * numpy.exp(3 * u_one * nu * x_coords[x_coords < 0]/(v*v))
    # Downstream
    outputArray[x_coords > 0] = N

    return outputArray
    """

    plot_properties_x = plot_properties.replace(
        interpol=1, legend_location="TopLeft"
    )
    plot_over_line_x, layout_x, line_chart_view_x = vfp.plot_f_lms_over_x(
        f_000_ana_x,
        results_folder,
        "spatial-distribution",
        plot_properties_x,
        direction="x",
        offset=[0.0, np.log(p_hat), 0.0],
        x_range=[bounds[0] * 0.1, bounds[1]],  # 10% of the upstream
    )
    # endregion

    return locals()


if __name__ in ["__main__", "__vtkconsole__"]:
    results = main()
    # Make all results available as global variables in a vtkconsole
    globals().update(results)
