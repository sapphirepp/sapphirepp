"""Plot the oblique shock example."""

import numpy as np
import matplotlib.pyplot as plt
import paraview.servermanager
from sapphireppplot import vfp, numpyify

plt.rcParams["figure.dpi"] = 150
plt.rcParams["figure.constrained_layout.use"] = True
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["figure.figsize"] = [8, 3]


def main() -> dict:
    """Plot the oblique shock example."""
    # ParaView setup and extraction of solution data
    plot_properties = vfp.PlotPropertiesVFP(
        dimension=2,
        momentum=True,
        scaled_distribution_function=True,
        debug_input_functions=True,
        sampling_resolution=1e-8,
    )

    results_folder, prm, solution, animation_scene = vfp.load_solution(
        plot_properties,
        path_prefix="$SAPPHIREPP_RESULTS/steady-state-oblique-shock",
    )

    # ParaView Plots
    # Spectrum at the shock
    (
        plot_over_line_p_at_the_shock,
        layout_p_at_the_shock,
        line_chart_view_p_at_the_shock,
    ) = vfp.plot_f_lms_over_p(
        solution,
        results_folder,
        "spectrum-at-shock",
        plot_properties,
        lms_indices=[[0, 0, 0]],
        offset=[0, 0, 0],
        value_range=[1e-5, 20],
    )

    # Spectrum far downstream
    downstream_position = 8000
    (
        plot_over_line_p_downstream,
        layout_p_downstream,
        line_chart_view_p_downstream,
    ) = vfp.plot_f_lms_over_p(
        solution,
        results_folder,
        "spectrum-downstream",
        plot_properties,
        lms_indices=[[0, 0, 0]],
        offset=[downstream_position, 0, 0],
        value_range=[1e-5, 20],
    )

    # Python plots

    # Skip numpy and matplotlib part for remote sessions
    if paraview.servermanager.ActiveConnection.IsRemote():
        return locals()

    def p_crit(diffusion_length: float, eta: float) -> float:
        u_sh = float(prm["Physical parameters"]["u_sh"])  # shock velocity
        magnetic_field_0 = float(
            prm["Physical parameters"]["B0"]
        )  # upstream magnetic field
        obliqueness = (
            float(prm["Physical parameters"]["obliqueness"]) / 180.0 * np.pi
        )
        return (
            diffusion_length
            * (3 * eta * u_sh * magnetic_field_0 * (eta**2 + 1))
            / (eta**2 + np.cos(obliqueness) ** 2)
        )

    p_inj = float(prm["Physical parameters"]["p_inj"])

    ln_p, g_000_at_the_shock = numpyify.to_numpy_1d(
        plot_over_line_p_at_the_shock,
        array_names=["g_000"],
        x_direction=1,
        x_min=np.log(p_inj),
    )

    ln_p, g_000_downstream = numpyify.to_numpy_1d(
        plot_over_line_p_downstream,
        array_names=["g_000"],
        x_direction=1,
        x_min=np.log(p_inj),
    )

    p_coordinate = np.exp(ln_p)

    # Iroshnikov-Kraichnan MHD turbulence: Spectra at the shock and far
    # downstream
    if (
        prm["Physical parameters"]["enhanced scattering zone"].lower()
        == "false"
    ):
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.loglog(
            p_coordinate,
            p_coordinate * g_000_at_the_shock[0, :],
            label="x = 0",
        )

        ax.loglog(
            p_coordinate,
            p_coordinate * g_000_downstream[0, :],
            label=f"x = {downstream_position:.2f}",
        )

        ax.legend(loc="upper left")
        ax.set_ylabel(r"$p^{4}f_{0}$")
        ax.set_xlabel(r"$p/mc$")
        ax.grid(True, which="both", alpha=0.5)
        plt.savefig(results_folder + "/spectra_kraichnan_turbulence.png")
    else:
        # Plots of spectra in the case of an enhanced scattering zone
        # Scattering regimes (a) and (b)
        case = prm["Physical parameters"]["case identifier"]
        x_min = -50
        x_max = 10
        p_value = 1

        plot_over_line_x, layout_x, line_chart_view_x = vfp.plot_f_lms_over_x(
            solution,
            results_folder,
            "spatial-profile",
            plot_properties,
            lms_indices=[[0, 0, 0], [1, 0, 0]],
            direction="x",
            offset=[0, np.log(p_value), 0],
            x_label=r"$x$",
            x_range=[x_min, x_max],
        )

        x, data = numpyify.to_numpy_1d(
            plot_over_line_x,
            array_names=["func_B_x", "func_B_z", "func_nu"],
            x_direction=0,
        )

        magnetic_field_magnitude = np.sqrt(data[0, :] ** 2 + data[1, :] ** 2)
        eta = data[2, :] / magnetic_field_magnitude * p_value
        fig, axs = plt.subplots(1, 2)

        axs[0].plot(
            x,
            data[2, :],
            label=r"$\nu(x, p = {p_value})$".format(p_value=p_value),
        )
        axs[0].legend(loc="upper left")
        axs[0].axvline(0, color="grey", linestyle="dashed")
        axs[0].set_xlabel(r"$x$")

        axs[1].plot(x, eta, label=r"$\eta(x)$")
        axs[1].legend(loc="upper left")
        axs[1].axvline(0, color="grey", linestyle="dashed")
        axs[1].set_xlabel(r"$x$")
        plt.savefig(results_folder + f"/scattering_regime_case_{case}.png")

        # The spectrum at the shock corresponding to either case (a) or case(b)

        # The scattering frequency inside the buffer (zeta) zone determines the
        # diffusion length, and thus the critical momentum for which the
        # diffusion length equals the size of the buffer zone
        p_crit_for_buffer_zone = p_crit(
            np.abs(float(prm["Physical parameters"]["transition point"])),
            float(prm["Physical parameters"]["zeta two"]),
        )

        fig, ax = plt.subplots()

        ax.grid(True, which="both", alpha=0.5)
        ax.set_xlabel(r"$p/mc$")
        ax.set_ylabel(r"$p^{4} f_{0}$")

        # Spectrum
        ax.loglog(
            p_coordinate,
            p_coordinate * g_000_at_the_shock[0, :],
            label=(
                r"$\eta = \zeta$"
                if case == "a"
                else r"$\nu = \omega_{\mathrm{g, up}} \zeta$"
            ),
        )

        # p_crit
        ax.axvline(p_crit_for_buffer_zone, color="grey", linestyle="dotted")
        ax.annotate(
            r"$p_{\text{crit}}$",
            xy=(p_crit_for_buffer_zone, ax.get_ylim()[0]),
            xytext=(
                (1.0 - 0.05) * p_crit_for_buffer_zone,
                (1.0 - 0.05) * ax.get_ylim()[0],
            ),
        )

        ax.legend()
        ax.set_xlim(left=1.5 * p_inj)
        plt.savefig(
            results_folder + f"/spectrum_enhanced_scattering_case_{case}.png"
        )

    return locals()


if __name__ in ["__main__", "__vtkconsole__"]:
    results = main()
    # Make all results available as global variables in a vtkconsole
    globals().update(results)
