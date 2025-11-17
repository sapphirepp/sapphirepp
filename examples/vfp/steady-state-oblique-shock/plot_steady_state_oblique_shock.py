"""Plot the oblique shock example."""

import numpy as np
from sapphireppplot import vfp, numpyify
import matplotlib.pyplot as plt

plt.rcParams["figure.dpi"] = 150
plt.rcParams["figure.constrained_layout.use"] = True
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["figure.figsize"] = [8, 3]


def main() -> dict:
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

    def p_crit(L_diff, eta):
        u_sh = float(prm["Physical parameters"]["u_sh"])  # shock velocity
        B = float(prm["Physical parameters"]["B0"])  # upstream magnetic field
        obliqueness = (
            float(prm["Physical parameters"]["obliqueness"]) / 180.0 * np.pi
        )
        return (
            L_diff
            * (3 * eta * u_sh * B * (eta**2 + 1))
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
            linewidth=2,
        )

        ax.loglog(
            p_coordinate,
            p_coordinate * g_000_downstream[0, :],
            label="x = {position:.2f}".format(position=downstream_position),
            linewidth=2,
        )

        ax.legend(loc="upper left")
        ax.set_ylabel(r"$p^{4}f_{0}$")
        ax.set_xlabel(r"$p/mc$")
        ax.grid(True, which="both", alpha=0.5)
        plt.savefig(results_folder + "/spectra_kraichnan_turbulence.png")
    else:
        # Plots of spectra in the case of an enhanced scattering zone
        # Scattering regimes (a) and (b)
        x_min = -50
        x_max = 10
        p_value = 1

        plot_over_line_x, layout_x, line_chart_view_x = vfp.plot_f_lms_over_x(
            solution,
            results_folder,
            "spatial-profile",
            plot_properties,
            lms_indices=[[0, 0, 0], [1, 0, 0]],
            direction=[
                [x_min, np.log(p_value), 0],
                [x_max, np.log(p_value), 0],
            ],
            x_label="arc length",
        )

        x, data = numpyify.to_numpy_1d(
            plot_over_line_x,
            array_names=["func_B_x", "func_B_z", "func_nu"],
            x_direction=0,
        )

        B_mag = np.sqrt(data[0, :] ** 2 + data[1, :] ** 2)
        eta = data[2, :] / B_mag * p_value
        fig, axs = plt.subplots(1, 2)

        axs[0].plot(x, data[2, :], label=r"$\nu(x, p = {p_value})$".format(p_value=p_value))
        axs[0].legend(loc="upper left")
        axs[0].axvline(0, color="grey", linestyle="dashed")
        axs[0].set_xlabel(r"$x$")

        axs[1].plot(x, eta, label=r"$\eta(x)$")
        axs[1].legend(loc="upper left")
        axs[1].axvline(0, color="grey", linestyle="dashed")
        axs[1].set_xlabel(r"$x$")
        plt.savefig(
            results_folder
            + "/scattering_regime_case_{case}.png".format(
                case=prm["Physical parameters"]["case identifier"]
            )
        )

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
                if prm["Physical parameters"]["case identifier"] == "a"
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
            results_folder
            + "/spectrum_enhanced_scattering_case_{case}.png".format(
                case=prm["Physical parameters"]["case identifier"]
            )
        )

    return locals()


if __name__ in ["__main__", "__vtkconsole__"]:
    results = main()
    # Make all results available as global variables in a vtkconsole
    globals().update(results)
