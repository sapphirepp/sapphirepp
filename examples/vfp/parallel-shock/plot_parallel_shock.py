"""Plot parallel-shock example."""

import numpy as np
import paraview.simple as ps
from sapphireppplot import vfp, pvplot, transform


def main() -> dict:
    """Plot parallel-shock example."""
    plot_properties = vfp.PlotPropertiesVFP(
        dimension=2,
        momentum=True,
        lms_indices=[[0, 0, 0], [1, 0, 0]],
        preview_size_2d=[1280, 720],
        camera_view_2d=(True, 0.85),
        color_bar_position=[0.1, 0.6],
        axes_stretch=[1.0, 16.0, 1.0],
    )

    results_folder, prm, solution, animation_scene = vfp.load_solution(
        plot_properties,
        path_prefix="$SAPPHIREPP_RESULTS/parallel-shock",
    )

    # region Plot 2D render view
    layout_2d, render_view_2d = vfp.plot_f_lms_2d(
        solution,
        results_folder,
        "parallel-shock-2D",
        plot_properties,
        value_range=[1e-6, 30.0],
        save_animation=True,
    )
    # endregion

    # region Plot shock region
    plot_properties_shock = plot_properties.replace(
        show_grid=True,
        axes_stretch=[1.0, 1.0, 1.0],
        text_color=[1.0, 1.0, 1.0],  # Changed color of label text of color bar
    )
    clip_shock = transform.clip_area(solution, -6.0, 6.0, plot_properties_shock)
    layout_shock_region, render_view_shock_region = vfp.plot_f_lms_2d(
        clip_shock,
        results_folder,
        "shock-region",
        plot_properties_shock,
        value_range=[1e-6, 30.0],
    )
    # endregion

    # region Physical parameters for the analytical solution
    source = float(prm["Physical parameters"]["Q"])  # type: ignore
    p_inj = float(prm["Physical parameters"]["p_inj"])  # type: ignore
    r = float(prm["Physical parameters"]["compression ratio"])  # type: ignore
    u_one = float(prm["Physical parameters"]["u_sh"])  # type: ignore
    nu_0 = float(prm["Physical parameters"]["nu0"])  # type: ignore
    delta_t = float(
        prm["VFP"]["Time stepping"]["Time step size"]
    )  # type: ignore
    p_hat = 10  # spatial dependence of f_000 is shown at p = p_hat
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
    outputArray = numpy.sqrt(4 * numpy.pi) * 3*source/(u_one * p_inj) \
        * r/(r - 1) * (p_coords/p_inj)**(-3*r/(r-1))

    return outputArray
    """

    plot_properties_p = plot_properties.replace(interpol=True)
    plot_over_line_p, layout_p, line_chart_view_p = vfp.plot_f_lms_over_p(
        f_000_ana_p,
        results_folder,
        "particle-spectrum",
        plot_properties_p,
        lms_indices=[[0, 0, 0]],
        offset=[0.001, 0, 0],
        x_range=[p_inj * -0.1, bounds[3]],
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
    # Normalization of the distribution function
    N = numpy.sqrt(4 * numpy.pi) * 3*source/(u_one * p_inj) \
        * r/(r - 1) * (p_hat/p_inj)**(-3*r/(r-1))

    x_coords = inputs[0].Points[:,0]
    outputArray = numpy.zeros(x_coords.shape[0])

    # Upstream
    outputArray[x_coords < 0] = N * numpy.exp(3 * u_one * nu_0 * x_coords[x_coords < 0]/(v*v))
    # Downstream
    outputArray[x_coords > 0] = N

    return outputArray
    """

    f_100_ana_x = ps.PythonCalculator(
        registrationName="f_100_ana_x", Input=f_000_ana_x
    )
    f_100_ana_x.ArrayName = "interpol_f_100"
    f_100_ana_x.UseMultilineExpression = 1
    f_100_ana_x.MultilineExpression = f"""
    # Physical Parameters
    source = {source}
    p_inj = {p_inj}
    r = {r}
    u_one = {u_one}
    nu_0 = {nu_0}
    p_hat = {p_hat}

    # Derived quantities
    v = p_hat/numpy.sqrt(p_hat**2 + 1)   # magnitude of the particle velocity
    # Normalization of the distribution function
    N = - numpy.sqrt(4 * numpy.pi/3) * 3 * u_one/v *  3*source/(u_one * p_inj) \
        * r/(r - 1) * (p_hat/p_inj)**(-3*r/(r-1))

    x_coords = inputs[0].Points[:,0]
    outputArray = numpy.zeros(x_coords.shape[0])

    # Upstream
    outputArray[x_coords < 0] = N * numpy.exp(3 * u_one * nu_0 * x_coords[x_coords < 0]/(v*v))
    # Downstream
    outputArray[x_coords > 0] = 0

    return outputArray
    """

    plot_properties_x = plot_properties.replace(
        interpol=1, legend_location="TopLeft"
    )
    plot_over_line_x, layout_x, line_chart_view_x = vfp.plot_f_lms_over_x(
        f_100_ana_x,
        results_folder,
        "spatial-distribution",
        plot_properties_x,
        lms_indices=[[0, 0, 0], [1, 0, 0]],
        direction="x",
        offset=[0.0, np.log(p_hat), 0.0],
    )
    # endregion

    # region Plot over time (time-dependent acceleration)
    probe_location, plot_properties_t = transform.probe_location(
        solution,
        [0.001, np.log(p_hat), 0.0],
        plot_properties_in=plot_properties.replace(
            legend_location="TopLeft", interpol=1
        ),
    )

    plot_over_time, plot_properties_t = transform.plot_over_time(
        probe_location,
        t_axes_scale=delta_t,
        results_folder=results_folder,
        filename="temporal-evolution",
        plot_properties_in=plot_properties_t,
    )

    f_000_ana_t = ps.PythonCalculator(
        registrationName="f_000_ana_t", Input=plot_over_time
    )
    f_000_ana_t.ArrayName = "interpol_f_000"
    f_000_ana_t.UseMultilineExpression = 1
    f_000_ana_t.MultilineExpression = f"""
    # Physical Parameters
    source = {source}
    p_inj = {p_inj}
    r = {r}
    u_one = {u_one}
    nu_0 = {nu_0}
    p_hat = {p_hat}
    # Normalization of the distribution function
    N = numpy.sqrt(4 * numpy.pi) * 3*source/(u_one * p_inj) \
        * r/(r - 1) * (p_hat/p_inj)**(-3*r/(r-1))

    # Work around to get a numpy array with the time steps
    # No idea if there is a better way
    t_coords = numpy.zeros((1,Time.GetSize()), dtype=np.float64)
    t_coords[:] = Time.GetArrays()
    t_coords = t_coords.flatten()[1:] # [1:] removes t = 0

    c1 = r / ( 2 * u_one**2 * nu_0) * (r + 1)/(r -1) \
        * numpy.log((1 + p_hat**2)/(1 + p_inj**2))
    c2 = 1/(3 * nu_0**2) * r/u_one**4 * (r**3 + 1)/(r - 1) \
        * (1/(1 + p_hat**2) - 1/(1 + p_inj**2)
            + numpy.log((1 + p_hat**2)/(1 + p_inj**2)))

    from scipy.special import erfc
    # phi_t is zero at t = 0
    phi_t = np.zeros(Time.GetSize())
    phi_t[1:] = 0.5 * (numpy.exp(2 * c1**2/c2) * erfc(numpy.sqrt(c1**3/(2 * t_coords * c2)) + numpy.sqrt(c1 * t_coords/(2*c2))) +erfc(numpy.sqrt(c1**3/(2 * t_coords * c2)) -  numpy.sqrt(c1 * t_coords/(2 * c2))))

    outputArray = N * phi_t

    return outputArray
    """

    layout_t = ps.CreateLayout("temporal-evolution")
    line_chart_view_t = pvplot.plot_line_chart_view(
        f_000_ana_t,
        layout_t,
        # x_label=r"$t$",
        x_label=r"$t / \Delta t$",
        y_label=r"$f_{lms}$",
        # x_array_name="Time",
        x_array_name="scaled_t_axes",
        visible_lines=["f_000 (id=0)", "f_100 (id=0)", "interpol_f_000 (id=0)"],
        log_y_scale=False,
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
