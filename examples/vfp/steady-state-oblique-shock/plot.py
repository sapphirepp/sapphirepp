import numpy as np
from sapphireppplot import vfp, numpyify, utils
import matplotlib.pyplot as plt

plot_properties = vfp.PlotPropertiesVFP(
    dimension=2,
    momentum=True,
    scaled_distribution_function=True,
    debug_input_functions=True,
    sampling_resolution=1e-8,
)

results_folder, prm, solution, animation_scene = vfp.load_solution(
    plot_properties,
    path_prefix="~/",
    results_folder="oblique-shock-example",
)

x_min = -50
x_max = 10
p_value = 20

plot_over_line_x, layout_x, line_chart_view_x = vfp.plot_f_lms_over_x(
    solution,
    results_folder,
    "spatial-profile",
    plot_properties,
    lms_indices=[[0, 0, 0], [1, 0, 0]],
    direction=[[x_min, np.log(p_value), 0], [x_max, np.log(p_value), 0]],
    x_label="arc length",
)

x, data = numpyify.to_numpy_1d(
    plot_over_line_x,
    array_names=["func_B_x", "func_B_z", "func_nu", "g_000"],
    x_direction=0,
)

B_mag = np.sqrt(data[0, :] ** 2 + data[1, :] ** 2)
eta = data[2, :] / B_mag * p_value
fig, axs = plt.subplots(1, 2, layout="constrained", figsize=(20, 10))

axs[0].plot(x, data[2, :], label=r"$\nu(x, p = 20)$")
axs[0].legend(loc="upper left", fontsize=12)
axs[0].axvline(0, color="grey", linestyle="dashed")
axs[0].set_xlabel(r"$x$", fontsize=12)

axs[1].plot(x, eta, label=r"$\eta(x)$")
axs[1].legend(loc="upper left", fontsize=12)
axs[1].axvline(0, color="grey", linestyle="dashed")
axs[1].set_xlabel(r"$x$", fontsize=12)


plot_over_line_p_at_the_shock, layout_p_at_the_shock, line_chart_view_p_at_the_shock = (
    vfp.plot_f_lms_over_p(
        solution,
        results_folder,
        "spectrum-at-shock",
        plot_properties,
        lms_indices=[[0, 0, 0]],
        offset=[0, 0, 0],
        value_range=[1e-5, 20],
    )
)

downstream_position = 50
plot_over_line_p_downstream, layout_p_downstream, line_chart_view_p_downstream = (
    vfp.plot_f_lms_over_p(
        solution,
        results_folder,
        "spectrum-downstream",
        plot_properties,
        lms_indices=[[0, 0, 0]],
        offset=[downstream_position, 0, 0],
        value_range=[1e-5, 20],
    )
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
fig, ax = plt.subplots(layout="constrained", figsize=(20, 10))

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

ax.legend(loc="upper right", fontsize=28)
ax.set_ylabel(r"$p^{4}f_{0}$", fontsize=32)
ax.set_xlabel(r"$p/mc$", fontsize=32)
# x_min = np.exp(np.fromstring(prm["VFP"]["Mesh"]["Point 1"], dtype=float, sep=",")[1])
# x_max = np.exp(np.fromstring(prm["VFP"]["Mesh"]["Point 2"], dtype=float, sep=",")[1])
# ax.set_ylim(, 300)
# ax.set_xlim(x_min, x_max)
# ax.set_xticks([1, 1e2, 1e4])
# ax.set_yticks([100, 300])
ax.tick_params(axis="both", which="both", labelsize=24)
ax.grid(True, which="both", alpha=0.5)

plt.show()
