import os
from sapphireppplot import vfp, pvload

results_folder = "/home/goyal/Sapphire_Code/sapphirepp/build/results/radiation-reaction-without-rotation"
figures_folder = os.path.join(results_folder, "figures_over_p")
os.makedirs(figures_folder, exist_ok=True)

lms_indices = [[0, 0, 0], [1, 0, 0]]

solution = pvload.load_solution_pvtu(
    results_folder,
    base_file_name="solution"
)

solution = pvload.scale_time_steps(solution, t_start=0.0, t_end=1.0)


for time_index, t in enumerate(solution.TimestepValues):

    solution.UpdatePipeline(t)

    plot_over_line, layout, line_chart_view = vfp.plot_f_lms_over_p(
        solution,
        figures_folder,
        f"f_vs_p_{time_index:04d}", 
        vfp.PlotPropertiesVFP(
            dimension=1,
            momentum=True
        ),
        lms_indices=lms_indices,
        direction="y",     
        x_range=[13.815510558, 18.420680744],  
        value_range=[1e-8, 1],          
        log_x_scale=False,           
        log_y_scale=True,            
        save_animation=True
    )

print("Done: PNGs saved in:", figures_folder)
