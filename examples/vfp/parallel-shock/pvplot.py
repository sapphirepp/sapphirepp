# trace generated using paraview version 5.13.2
# import paraview
# paraview.compatibility.major = 5
# paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


# =============================================================================
# Define results path
# =============================================================================

# Either using command line argument, or using input
results_folder = ""
if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        results_folder = sys.argv[1]
if not results_folder:
    results_folder = input("Path to results folder: ")
print(f"Using results in '{results_folder}'")


# =============================================================================
# Load pvtu files
# =============================================================================

from paraview.util import *

pvtu_files = paraview.util.Glob(results_folder + "/solution_*.pvtu")
if not pvtu_files:
    raise FileNotFoundError(
        f"No .pvtu files found matching '{results_folder}/solution_*.pvtu'"
    )

# create a new 'XML Partitioned Unstructured Grid Reader'
solution = XMLPartitionedUnstructuredGridReader(
    registrationName="solution",
    FileName=pvtu_files,
)
# for pvtu, the time variable does not work
solution.TimeArray = "None"
# Only show 'f_000' and 'f_100'
solution.PointArrayStatus = ["f_000", "f_100"]


# # =============================================================================
# # Plot 2D render view
# # =============================================================================

# # get animation scene
# animationScene1 = GetAnimationScene()

# # update animation scene based on data timesteps
# animationScene1.UpdateAnimationUsingDataTimeSteps()

# # create new layout object '2D Plot'
# # Using a new layout for every plot simplifies the scaling of fonts
# layout2D = CreateLayout(name="2D Plot")

# # Create a new 'Render View'
# renderView1 = CreateView("RenderView")

# # add view to a layout so it's visible in UI
# AssignViewToLayout(view=renderView1, layout=layout2D, hint=0)

# # set active view
# SetActiveView(renderView1)

# # show data in view
# solutionDisplay = Show(solution, renderView1, "UnstructuredGridRepresentation")

# # PreviewMode and SetSize must be done after Show
# # Enter preview mode
# layout2D.PreviewMode = [1024, 1024]
# # layout/tab size in pixels
# layout2D.SetSize(1024, 1024)

# # reset view to fit data
# renderView1.ResetCamera(False, 0.9)

# # Hide orientation axes
# renderView1.OrientationAxesVisibility = 0

# # current camera placement for renderView1
# renderView1.InteractionMode = "2D"

# # Properties modified on renderView1
# renderView1.UseColorPaletteForBackground = 0
# renderView1.BackgroundColorMode = "Single Color"
# renderView1.Background = [1.0, 1.0, 1.0]

# # Properties modified on solutionDisplay
# solutionDisplay.DisableLighting = 1
# solutionDisplay.Diffuse = 1.0


# # ------------------
# # Show DataAxesGrid
# # ------------------

# solutionDisplay.DataAxesGrid.GridAxesVisibility = 1
# solutionDisplay.DataAxesGrid.XTitle = "$x$"
# solutionDisplay.DataAxesGrid.YTitle = "$\\ln p$"
# # Only show Axes Min-X//Y/Z
# solutionDisplay.DataAxesGrid.AxesToLabel = 7
# # solutionDisplay.DataAxesGrid.FacesToRender = 7
# # Set default font size: 24 for title, 18 for label
# solutionDisplay.DataAxesGrid.XTitleFontSize = 24
# solutionDisplay.DataAxesGrid.YTitleFontSize = 24
# solutionDisplay.DataAxesGrid.XLabelFontSize = 18
# solutionDisplay.DataAxesGrid.YLabelFontSize = 18
# # Use gray color for label for good visibility in both light and dark mode
# solutionDisplay.DataAxesGrid.XTitleColor = [0.5, 0.5, 0.5]
# solutionDisplay.DataAxesGrid.YTitleColor = [0.5, 0.5, 0.5]
# solutionDisplay.DataAxesGrid.XLabelColor = [0.5, 0.5, 0.5]
# solutionDisplay.DataAxesGrid.YLabelColor = [0.5, 0.5, 0.5]
# solutionDisplay.DataAxesGrid.GridColor = [0.5, 0.5, 0.5]


# # update the view to ensure updated data information
# renderView1.Update()


# # ----------------------
# # Create color bar plot
# # ----------------------

# # set scalar coloring
# ColorBy(solutionDisplay, ("POINTS", "f_000"))

# # show color bar/color legend
# solutionDisplay.SetScalarBarVisibility(renderView1, True)

# # get color transfer function/color map for 'f_000'
# f_000LUT = GetColorTransferFunction("f_000")
# # get opacity transfer function/opacity map for 'f_000'
# f_000PWF = GetOpacityTransferFunction("f_000")
# # get 2D transfer function for 'f_000'
# f_000TF2D = GetTransferFunction2D("f_000")

# # Rescale transfer function
# f_000LUT.RescaleTransferFunction(1e-2, 10.0)
# # Rescale transfer function
# f_000PWF.RescaleTransferFunction(1e-2, 10.0)
# # Rescale 2D transfer function
# f_000TF2D.RescaleTransferFunction(1e-2, 10.0, 0.0, 1.0)

# # convert to log space
# f_000LUT.MapControlPointsToLogSpace()

# # Properties modified on f_000LUT
# f_000LUT.UseLogScale = 1

# # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
# f_000LUT.ApplyPreset("Viridis (matplotlib)", True)

# # get color legend/bar for f_000LUT in view renderView1
# f_000LUTColorBar = GetScalarBar(f_000LUT, renderView1)

# # Properties modified on f_000LUTColorBar
# f_000LUTColorBar.Title = "$f_{000}$"
# f_000LUTColorBar.TitleFontSize = 24
# f_000LUTColorBar.LabelFontSize = 18
# f_000LUTColorBar.TitleColor = [0.5, 0.5, 0.5]
# f_000LUTColorBar.LabelColor = [0.5, 0.5, 0.5]

# # change scalar bar placement
# f_000LUTColorBar.WindowLocation = "Any Location"
# f_000LUTColorBar.ScalarBarLength = 0.25
# f_000LUTColorBar.Position = [0.15, 0.55]


# # Go to last simulation time step
# animationScene1.GoToLast()


# # ----------------
# # Save screenshot
# # ----------------

# # Save screenshot in results folder (at data server).
# print(f"Save screenshot '{results_folder}/quick-start-2D.png'")

# # save screenshot
# SaveScreenshot(
#     filename=results_folder + "/parallel-shock-2D.png",
#     viewOrLayout=renderView1,
#     location=vtkPVSession.DATA_SERVER,
#     TransparentBackground=1,
# )


# # ----------------
# # Save Animation
# # ----------------

# print(f"Save animation '{results_folder}/parallel-shock-2D.*.png'")

# # save animation
# SaveAnimation(
#     filename=results_folder + "/parallel-shock-2D.png",
#     viewOrLayout=renderView1,
#     location=vtkPVSession.DATA_SERVER,
#     TransparentBackground=1,
#     FrameStride=1,
# )

# =============================================================================
# Physical Parameters for the analytical solution
# =============================================================================
Q = 1.              # number density of the injected particles
p_inj = 2.           # magnitude of the injection momentum
r = 4                # compression ratio of the shock
u_one = 0.0167       # shock velocity
nu_0 = 1             # inverse of the hall parameter h, i.e. h = omega_g/nu
p_hat = 50           # spatial dependence of f_000 is shown at p = p_hat

# =============================================================================
# Plot f(x)
# =============================================================================
final_solution = ExtractTimeSteps(solution)

width = 1280
height = 720
# Get the bounds in x
# Fetch data information from the solution
solution_data = servermanager.Fetch(solution)
# Get bounds of the data
bounds = solution_data.GetBounds()

# -----------------------
# Compute the analytical solution
# ----------------------
# final_solution = ExtractTimeSteps(solution, )
# python_calc_ana_x = PythonCalculator(registrationName="f_000_ana_x", Input=solution)
# python_calc_ana_x.ArrayName = 'f_000_ana_x'
# python_calc_ana_x.UseMultilineExpression = 1

# python_cal_expression_x = """#########
# # Physical Parameters
# Q = {Q}              # number density of the injected particles
# p_inj = {p_inj}      # magnitude of the injection momentum
# r = {r}              # compression ratio of the shock
# u_one = {u_one}      # shock velocity
# nu_0 = {nu_0}        # inverse of the hall parameter h, i.e. h = omega_g/nu
# p_hat = {p_hat}      # spatial dependence of f_000 at p_hat

# # Derived quantities
# v = p_hat/numpy.sqrt(p_hat**2 + 1)   # magnitude of the particle velocity
# nu = nu_0/p_hat                      # scattering frequency for p_hat
# # Normalization of the distribution function
# N = 3*Q/(numpy.sqrt(4 * numpy.pi) * u_one * p_inj**3) \
#     * r/(r - 1) * (p_hat/p_inj)**(-3*r/(r-1))  

# x_coords = inputs[0].Points[:,0]
# outputArray = numpy.zeros(x_coords.shape[0])

# # Upstream
# outputArray[x_coords < 0] = N * numpy.exp(3 * u_one * nu * x_coords[x_coords < 0]/(v*v))
# # Downstream
# outputArray[x_coords > 0] = N

# return outputArray""".format(Q = Q,
#                              p_inj = p_inj,
#                              r = r,
#                              u_one = u_one,
#                              nu_0 = nu_0,
#                              p_hat = p_hat)

# python_calc_ana_x.MultilineExpression = python_cal_expression_x

# -----------------------
# Create a 'Plot Over Line'
# -----------------------

plotOverLine_f_x = PlotOverLine(registrationName="f(x) Plot", Input=solution)

# Extract min_x and max_x from the bounds
min_x = bounds[0] * 0.1        # 10 % of the upstream
max_x = bounds[1]

from math import log
ln_p_hat = log(p_hat)

# Specify where to plot over line
plotOverLine_f_x.Point1 = [min_x,ln_p_hat, 0.0]
plotOverLine_f_x.Point2 = [max_x,ln_p_hat, 0.0]

# Set SamplingPattern
plotOverLine_f_x.SamplingPattern = "Sample At Segment Centers"
# The automatic ComputeTolerance does not work for large shock-grids,
# so we set it manually
plotOverLine_f_x.ComputeTolerance = 0
# A tolerance of 1e-8 works well,
# this is approx (dx_min / dx_max) / 1000,
# with dx_min/max the minimum/maximum cell size in x
plotOverLine_f_x.Tolerance = 1e-8

# ------------------------------------
# Create new layout and LineChartView
# ------------------------------------

# create new layout object 'f(x) Plot'
layout_f_x = CreateLayout(name="f(x) Plot")

# Enter preview mode
layout_f_x.PreviewMode = [width, height]
# layout/tab size in pixels
layout_f_x.SetSize(width, height)

# Create a new 'Line Chart View'
lineChartView_f_x = CreateView("XYChartView")
# NOTE: ViewSize needs to be set if pvbatch(pvpython) is used to produce the
# images
lineChartView_f_x.ViewSize = [width,height] 
# lineChartView_f_x.ChartTitle = "Spatial dependence"
lineChartView_f_x.LegendLocation = 'TopLeft'
lineChartView_f_x.LeftAxisTitle = "$f_{000}(x)$"
lineChartView_f_x.LeftAxisLogScale = 0
lineChartView_f_x.BottomAxisTitle = "$x$"
lineChartView_f_x.ChartTitleFontSize = 30
lineChartView_f_x.LeftAxisTitleFontSize = 24
lineChartView_f_x.BottomAxisTitleFontSize = 24
lineChartView_f_x.LegendFontSize = 18
lineChartView_f_x.LeftAxisLabelFontSize = 18
lineChartView_f_x.BottomAxisLabelFontSize = 18

# assign view to a particular cell in the layout
AssignViewToLayout(view=lineChartView_f_x, layout=layout_f_x, hint=0)

# set active view
SetActiveView(lineChartView_f_x)

# set active source
SetActiveSource(plotOverLine_f_x)

# ---------------------
# Display PlotOverLine
# ---------------------

# show data in view
fxPlotDisplay = Show(
    plotOverLine_f_x, lineChartView_f_x, "XYChartRepresentation"
)

# Set which data to use for the x-axis
fxPlotDisplay.XArrayName = "Points_X"

# Set the labels of the plots
fxPlotDisplay.SeriesLabel =  [ 'f_000', '$f_{000}$', 'f_000_ana_x', '$f_{000}$-ana']

# Adapt line thickness
fxPlotDisplay.SeriesLineThickness = ['f_000', '2', 'f_000_ana_x', '2']

# Line style
fxPlotDisplay.SeriesLineStyle = ['f_000', '1', 'f_000_ana_x', '2'] # ana dashed

# Color the plots
fxPlotDisplay.SeriesColor = [
    'f_000', '0.10980392156862745', '0.5843137254901961', '0.8039215686274',
    'f_000_ana_x', '0.25882352941176473', '0.23921568627450981', '0.6627450980392157'
]

# Ensure that f_000 and f_000_ana_x are displayed
fxPlotDisplay.SeriesVisibility = ['f_000', 'f_000_ana_x']

# ----------------
# Save screenshot
# ----------------

print(f"Save screenshot '{results_folder}/spatial-distribution.png'")

# save screenshot
SaveScreenshot(
    filename=results_folder + '/spatial-distribution.png',
    location = vtkPVSession.DATA_SERVER,
    viewOrLayout=layout_f_x,
)


# ----------
# Save data
# ----------

print(f"Save data '{results_folder}/spatial-distribution.csv'")

# save data
SaveData(
    filename=results_folder + "/spatial-distribution.csv",
    proxy=plotOverLine_f_x,
    location = vtkPVSession.DATA_SERVER,
    ChooseArraysToWrite=2,
    PointDataArrays=['f_000','f_000_ana_x'],
    Precision=6,
    UseScientificNotation=1,
)
