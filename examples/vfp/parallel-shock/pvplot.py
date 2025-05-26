# trace generated using paraview version 5.13.2
# import paraview
# paraview.compatibility.major = 5
# paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *

#### disable automatic camera reset on 'Show'
# paraview.simple._DisableFirstRenderCameraReset()

# =============================================================================
# Define results path
# =============================================================================
results_folder = ""
# Either using command line argument, or using input
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
solution.UpdatePipelineInformation()

# for pvtu, the time variable does not work
# solution.TimeArray = "None"
# Only show 'f_000' and 'f_100'
solution.PointArrayStatus = ["f_000", "f_100"]

# =============================================================================
# Properties of the images of the plots
# =============================================================================
width = 1280
height = 720

# =============================================================================
# Plot 2D render view
# =============================================================================

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create new layout object '2D Plot'
# Using a new layout for every plot simplifies the scaling of fonts
layout2D = CreateLayout(name="2D Plot")

# Create a new 'Render View'
renderView1 = CreateView("RenderView")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=renderView1, layout=layout2D, hint=0)

# set active view
SetActiveView(renderView1)

# show data in view
solutionDisplay = Show(solution, renderView1, "UnstructuredGridRepresentation")

# PreviewMode and SetSize must be done after Show
# Enter preview mode
layout2D.PreviewMode = [width, height]
# layout/tab size in pixels
layout2D.SetSize(width, height)

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# current camera placement for renderView1
renderView1.InteractionMode = "2D"

# Properties modified on renderView1
renderView1.UseColorPaletteForBackground = 0
renderView1.BackgroundColorMode = "Single Color"
renderView1.Background = [1.0, 1.0, 1.0]

# Properties modified on solutionDisplay
solutionDisplay.DisableLighting = 1
solutionDisplay.Diffuse = 1.0

# --------------
# Show AxesGrid
# --------------

# We have to use AxesGrid instead of DataAxesGrid
# to support scaling of the p axes (see below)

renderView1.AxesGrid.Visibility = 1
renderView1.AxesGrid.XTitle = "$x$"
renderView1.AxesGrid.YTitle = "$\\ln p$"
# Only show Axes Min-X//Y/Z
renderView1.AxesGrid.AxesToLabel = 7
# renderView1.AxesGrid.FacesToRender = 7
# Set default font size: 24 for title, 18 for label
renderView1.AxesGrid.XTitleFontSize = 24
renderView1.AxesGrid.YTitleFontSize = 24
renderView1.AxesGrid.XLabelFontSize = 18
renderView1.AxesGrid.YLabelFontSize = 18
# Use gray color for label for good visibility in both light and dark mode
renderView1.AxesGrid.XTitleColor = [0.5, 0.5, 0.5]
renderView1.AxesGrid.YTitleColor = [0.5, 0.5, 0.5]
renderView1.AxesGrid.XLabelColor = [0.5, 0.5, 0.5]
renderView1.AxesGrid.YLabelColor = [0.5, 0.5, 0.5]
renderView1.AxesGrid.GridColor = [0.5, 0.5, 0.5]


# -------------
# Scale p axes
# -------------

# Define factor to scale p axes
scale_p_axes = 16.0

# Properties modified on solutionDisplay
solutionDisplay.Scale = [1.0, scale_p_axes, 1.0]

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.DataScale = [1.0, scale_p_axes, 1.0]


# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data
renderView1.ResetCamera(True, 0.85)


# ----------------------
# Create color bar plot
# ----------------------

# change use separate color map
solutionDisplay.UseSeparateColorMap = 1

# set scalar coloring
ColorBy(solutionDisplay, ("POINTS", "f_000"), True)

# show color bar/color legend
solutionDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'f_000'
f_000LUT = GetColorTransferFunction("f_000", solutionDisplay, separate=True)
# get opacity transfer function/opacity map for 'f_000'
f_000PWF = GetOpacityTransferFunction("f_000", solutionDisplay, separate=True)
# get 2D transfer function for 'f_000'
f_000TF2D = GetTransferFunction2D("f_000", solutionDisplay, separate=True)

# Rescale transfer function
f_000LUT.RescaleTransferFunction(1e-6, 30.0)
# Rescale transfer function
f_000PWF.RescaleTransferFunction(1e-6, 30.0)
# Rescale 2D transfer function
f_000TF2D.RescaleTransferFunction(1e-6, 30.0, 0.0, 1.0)

# convert to log space
f_000LUT.MapControlPointsToLogSpace()

# Properties modified on f_000LUT
f_000LUT.UseLogScale = 1

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
f_000LUT.ApplyPreset("Viridis (matplotlib)", True)

# get color legend/bar for f_000LUT in view renderView1
f_000LUTColorBar = GetScalarBar(f_000LUT, renderView1)

# Properties modified on f_000LUTColorBar
f_000LUTColorBar.Title = "$f_{000}$"
f_000LUTColorBar.TitleFontSize = 24
f_000LUTColorBar.LabelFontSize = 18
f_000LUTColorBar.TitleColor = [0.5, 0.5, 0.5]
f_000LUTColorBar.LabelColor = [0.5, 0.5, 0.5]

# change scalar bar placement
f_000LUTColorBar.WindowLocation = "Any Location"
f_000LUTColorBar.ScalarBarLength = 0.25
f_000LUTColorBar.Position = [0.1, 0.6]

# Go to last simulation time step
animationScene1.GoToLast()

# ----------------
# Save screenshot
# ----------------

# Save screenshot in results folder (at data server).
print(f"Save screenshot '{results_folder}/parallel-shock-2D.png'")

# save screenshot
SaveScreenshot(
    filename=results_folder + "/parallel-shock-2D.png",
    viewOrLayout=renderView1,
    location=vtkPVSession.DATA_SERVER,
    TransparentBackground=1,
)

# ----------------
# Save Animation
# ----------------

print(f"Save animation '{results_folder}/parallel-shock-2D.*.png'")

# save animation
SaveAnimation(
    filename=results_folder + "/parallel-shock-2D.png",
    viewOrLayout=renderView1,
    location=vtkPVSession.DATA_SERVER,
    TransparentBackground=1,
    FrameStride=1,
)


# =============================================================================
# Create Clip for shock region
# =============================================================================

# create a new 'Clip'
clipShock = Clip(registrationName="Clip Shock", Input=solution)

# toggle interactive widget visibility (only when running from the GUI)
ShowInteractiveWidgets(proxy=clipShock.ClipType)

# Properties modified on clipShock
clipShock.ClipType = "Box"
clipShock.Crinkleclip = 1

# Get the bounds of solution
# Fetch data information from the solution
solution_data = servermanager.Fetch(solution)
# Get bounds of the data
bounds = solution_data.GetBounds()
# Extract min_y and max_y from the bounds
min_y = bounds[2]
max_y = bounds[3]

# Use small epsilon in z to capture cells inside box
epsilon_z = 0.1

# Define min_x and max_x to capture shock
min_x = -6.0
max_x = 6.0

# Properties modified on clipShock.ClipType
clipShock.ClipType.Position = [min_x, min_y, -epsilon_z]
clipShock.ClipType.Length = [max_x - min_x, max_y - min_y, 2 * epsilon_z]

# toggle interactive widget visibility (only when running from the GUI)
HideInteractiveWidgets(proxy=clipShock.ClipType)


# =============================================================================
# Plot 2D render view of shock region
# =============================================================================

# create new layout object 'Shock Region'
layoutShock = CreateLayout(name="Shock Region")

# Create a new 'Render View'
renderViewShock = CreateView("RenderView")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=renderViewShock, layout=layoutShock, hint=0)

# set active view
SetActiveView(renderViewShock)

# show data in view
clipDisplay = Show(clipShock, renderViewShock, "UnstructuredGridRepresentation")

# change representation type
clipDisplay.SetRepresentationType("Surface With Edges")

# PreviewMode and SetSize must be done after Show
# Enter preview mode
layoutShock.PreviewMode = [width, height]
# layout/tab size in pixels
layoutShock.SetSize(width, height)

# reset view to fit data
renderViewShock.ResetCamera(False, 0.9)

# Hide orientation axes
renderViewShock.OrientationAxesVisibility = 0

# current camera placement for renderViewShock
renderViewShock.InteractionMode = "2D"

# Properties modified on renderViewShock
renderViewShock.UseColorPaletteForBackground = 0
renderViewShock.BackgroundColorMode = "Single Color"
renderViewShock.Background = [1.0, 1.0, 1.0]

# Properties modified on clipDisplay
clipDisplay.DisableLighting = 1
clipDisplay.Diffuse = 1.0


# --------------
# Show AxesGrid
# --------------

# We have to use AxesGrid instead of DataAxesGrid
# to support scaling of the p axes (see below)

renderViewShock.AxesGrid.Visibility = 1
renderViewShock.AxesGrid.XTitle = "$x$"
renderViewShock.AxesGrid.YTitle = "$\\ln p$"
# Only show Axes Min-X//Y/Z
renderViewShock.AxesGrid.AxesToLabel = 7
# renderViewShock.AxesGrid.FacesToRender = 7
# Set default font size: 24 for title, 18 for label
renderViewShock.AxesGrid.XTitleFontSize = 24
renderViewShock.AxesGrid.YTitleFontSize = 24
renderViewShock.AxesGrid.XLabelFontSize = 18
renderViewShock.AxesGrid.YLabelFontSize = 18
# Use gray color for label for good visibility in both light and dark mode
renderViewShock.AxesGrid.XTitleColor = [0.5, 0.5, 0.5]
renderViewShock.AxesGrid.YTitleColor = [0.5, 0.5, 0.5]
renderViewShock.AxesGrid.XLabelColor = [0.5, 0.5, 0.5]
renderViewShock.AxesGrid.YLabelColor = [0.5, 0.5, 0.5]
renderViewShock.AxesGrid.GridColor = [0.5, 0.5, 0.5]


# -------------
# Scale p axes
# -------------

# Define factor to scale p axes
scale_p_axes = 1.0
# Properties modified on clipDisplay
clipDisplay.Scale = [1.0, scale_p_axes, 1.0]

# Properties modified on renderViewShock.AxesGrid
renderViewShock.AxesGrid.DataScale = [1.0, scale_p_axes, 1.0]


# update the view to ensure updated data information
renderViewShock.Update()

# reset view to fit data
renderViewShock.ResetCamera(True, 0.85)


# ----------------------
# Create color bar plot
# ----------------------

# change use separate color map
clipDisplay.UseSeparateColorMap = 1

# set scalar coloring
ColorBy(clipDisplay, ("POINTS", "f_000"), True)

# show color bar/color legend
clipDisplay.SetScalarBarVisibility(renderViewShock, True)

# get color transfer function/color map for 'f_000'
f_000LUT = GetColorTransferFunction("f_000", clipDisplay, separate=True)
# get opacity transfer function/opacity map for 'f_000'
f_000PWF = GetOpacityTransferFunction("f_000", clipDisplay, separate=True)
# get 2D transfer function for 'f_000'
f_000TF2D = GetTransferFunction2D("f_000", clipDisplay, separate=True)

# Rescale transfer function
f_000LUT.RescaleTransferFunction(1e-6, 30.0)
# Rescale transfer function
f_000PWF.RescaleTransferFunction(1e-6, 30.0)
# Rescale 2D transfer function
f_000TF2D.RescaleTransferFunction(1e-6, 30.0, 0.0, 1.0)

# convert to log space
f_000LUT.MapControlPointsToLogSpace()

# Properties modified on f_000LUT
f_000LUT.UseLogScale = 1

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
f_000LUT.ApplyPreset("Viridis (matplotlib)", True)

# get color legend/bar for f_000LUT in view renderViewShock
f_000LUTColorBar = GetScalarBar(f_000LUT, renderViewShock)

# Properties modified on f_000LUTColorBar
f_000LUTColorBar.Title = "$f_{000}$"
f_000LUTColorBar.TitleFontSize = 24
f_000LUTColorBar.LabelFontSize = 18
f_000LUTColorBar.TitleColor = [0.8, 0.8, 0.8]
f_000LUTColorBar.LabelColor = [0.8, 0.8, 0.8]

# change scalar bar placement
f_000LUTColorBar.WindowLocation = "Any Location"
f_000LUTColorBar.ScalarBarLength = 0.25
f_000LUTColorBar.Position = [0.1, 0.6]


# ----------------
# Save screenshot
# ----------------

# Save screenshot in results folder (at data server).
print(f"Save screenshot '{results_folder}/shock-region.png'")

# save screenshot
SaveScreenshot(
    filename=results_folder + "/shock-region.png",
    viewOrLayout=renderViewShock,
    location=vtkPVSession.DATA_SERVER,
    TransparentBackground=1,
)


# =============================================================================
# Physical Parameters for the analytical solution
# =============================================================================
Q = 0.1  # number density of the injected particles
p_inj = 1.0  # magnitude of the injection momentum
r = 4  # compression ratio of the shock
u_one = 0.03  # shock velocity
nu_0 = 1  # inverse of the hall parameter h, i.e. h = omega_g/nu
p_hat = 10  # spatial dependence of f_000 is shown at p = p_hat
delta_t = 100

# =============================================================================
# Plots Over Line
# =============================================================================

# Get bounds of the data
bounds = solution.GetDataInformation().GetBounds()

# -----------------------
# Compute the analytical solution
# ----------------------
# =============================================================================
# Plot f(p)
# =============================================================================
python_calc_f_000_ana_p = PythonCalculator(
    registrationName="f_000_ana_p", Input=solution
)
python_calc_f_000_ana_p.ArrayName = "f_000_ana_p"
python_calc_f_000_ana_p.UseMultilineExpression = 1

python_calc_expression_f_000_p = """#########
# Physical Parameters
Q = {Q}              # number density of the injected particles
p_inj = {p_inj}      # magnitude of the injection momentum
r = {r}              # compression ratio of the shock
u_one = {u_one}      # shock velocity

p_coords = numpy.exp(inputs[0].Points[:,1]) 
outputArray = numpy.sqrt(4 * numpy.pi) * 3*Q/(u_one * p_inj) \
    * r/(r - 1) * (p_coords/p_inj)**(-3*r/(r-1))
return outputArray""".format(
    Q=Q, p_inj=p_inj, r=r, u_one=u_one
)

python_calc_f_000_ana_p.MultilineExpression = python_calc_expression_f_000_p

plotOverLine_f_p = PlotOverLine(
    registrationName="f(p) Plot", Input=python_calc_f_000_ana_p
)

# Set SamplingPattern
plotOverLine_f_p.SamplingPattern = "Sample At Segment Centers"

# Extract min_x and max_x from the bounds
min_y = p_inj * -0.1
max_y = bounds[3]

# Specify where to plot over line
plotOverLine_f_p.Point1 = [0.001, min_y, 0.0]
plotOverLine_f_p.Point2 = [0.001, max_y, 0.0]

# ------------------------------------
# Create new layout and LineChartView
# ------------------------------------

# create new layout object 'f(x) Plot'
layout_f_p = CreateLayout(name="f(p) Plot")

# Enter preview mode
layout_f_p.PreviewMode = [width, height]
# layout/tab size in pixels
layout_f_p.SetSize(width, height)

# Create a new 'Line Chart View'
lineChartView_f_p = CreateView("XYChartView")
# NOTE: ViewSize needs to be set if pvbatch(pvpython) is used to produce the
# images
lineChartView_f_p.ViewSize = [width, height]
lineChartView_f_p.ViewTime = solution.TimestepValues[-1]  # last time step
lineChartView_f_p.LegendLocation = "TopRight"
lineChartView_f_p.LeftAxisTitle = "$f_{000}(p)$"
lineChartView_f_p.LeftAxisLogScale = 1
lineChartView_f_p.BottomAxisTitle = "$\\ln p$"
lineChartView_f_p.ChartTitleFontSize = 30
lineChartView_f_p.LeftAxisTitleFontSize = 24
lineChartView_f_p.BottomAxisTitleFontSize = 24
lineChartView_f_p.LegendFontSize = 18
lineChartView_f_p.LeftAxisLabelFontSize = 18
lineChartView_f_p.BottomAxisLabelFontSize = 18

# assign view to a particular cell in the layout
AssignViewToLayout(view=lineChartView_f_p, layout=layout_f_p, hint=0)

# ---------------------
# Display PlotOverLine
# ---------------------

# show data in view
fpPlotDisplay = Show(
    plotOverLine_f_p, lineChartView_f_p, "XYChartRepresentation"
)
# Set which data to use for the x-axis
fpPlotDisplay.XArrayName = "Points_Y"

# Set the labels of the plots
fpPlotDisplay.SeriesLabel = [
    "f_000",
    "$f_{000}$",
    "f_000_ana_p",
    "$f_{000}$-ana",
]

# Adapt line thickness
fpPlotDisplay.SeriesLineThickness = ["f_000", "2", "f_000_ana_p", "2"]

# Line style
fpPlotDisplay.SeriesLineStyle = ["f_000", "1", "f_000_ana_p", "2"]  # ana dashed

# Color the plots
fpPlotDisplay.SeriesColor = [
    "f_000",
    "0.10980392156862745",
    "0.5843137254901961",
    "0.8039215686274",
    "f_000_ana_p",
    "0.25882352941176473",
    "0.23921568627450981",
    "0.6627450980392157",
]

fpPlotDisplay.SeriesVisibility = ["f_000", "f_000_ana_p"]

# ----------------
# Save screenshot
# ----------------

print(f"Save screenshot '{results_folder}/particle-spectrum.png'")

# save screenshot
SaveScreenshot(
    filename=results_folder + "/particle-spectrum.png",
    location=vtkPVSession.DATA_SERVER,
    viewOrLayout=layout_f_p,
)

# ----------
# Save data
# ----------

print(f"Save data '{results_folder}/particle-spectrum.csv'")

# save data
SaveData(
    filename=results_folder + "/particle-spectrum.csv",
    proxy=plotOverLine_f_p,
    location=vtkPVSession.DATA_SERVER,
    ChooseArraysToWrite=2,
    PointDataArrays=["f_000", "f_000_ana_p"],
    Precision=6,
    UseScientificNotation=1,
)

# =============================================================================
# Plot f(x)
# =============================================================================
# -----------------------
# Compute the analytical solution
# ----------------------

python_calc_f_000_ana_x = PythonCalculator(
    registrationName="f_000_ana_x", Input=solution
)
python_calc_f_000_ana_x.ArrayName = "f_000_ana_x"
python_calc_f_000_ana_x.UseMultilineExpression = 1

python_calc_expression_f_000_x = """#########
# Physical Parameters
Q = {Q}              # number density of the injected particles
p_inj = {p_inj}      # magnitude of the injection momentum
r = {r}              # compression ratio of the shock
u_one = {u_one}      # shock velocity
nu_0 = {nu_0}        # inverse of the hall parameter h, i.e. h = omega_g/nu
p_hat = {p_hat}      # spatial dependence of f_000 at p_hat

# Derived quantities
v = p_hat/numpy.sqrt(p_hat**2 + 1)   # magnitude of the particle velocity
# Normalization of the distribution functionx
N = numpy.sqrt(4 * numpy.pi) * 3*Q/(u_one * p_inj) \
    * r/(r - 1) * (p_hat/p_inj)**(-3*r/(r-1))

x_coords = inputs[0].Points[:,0]
outputArray = numpy.zeros(x_coords.shape[0])

# Upstream
outputArray[x_coords < 0] = N * numpy.exp(3 * u_one * nu_0 * x_coords[x_coords < 0]/(v*v))
# Downstream
outputArray[x_coords > 0] = N

return outputArray""".format(
    Q=Q, p_inj=p_inj, r=r, u_one=u_one, nu_0=nu_0, p_hat=p_hat
)

python_calc_f_000_ana_x.MultilineExpression = python_calc_expression_f_000_x

python_calc_f_100_ana_x = PythonCalculator(
    registrationName="f_100_ana_x", Input=python_calc_f_000_ana_x
)
python_calc_f_100_ana_x.ArrayName = "f_100_ana_x"
python_calc_f_100_ana_x.UseMultilineExpression = 1

python_calc_expression_f_100_x = """#########
# Physical Parameters
Q = {Q}              # number density of the injected particles
p_inj = {p_inj}      # magnitude of the injection momentum
r = {r}              # compression ratio of the shock
u_one = {u_one}      # shock velocity
nu_0 = {nu_0}        # inverse of the hall parameter h, i.e. h = omega_g/nu
p_hat = {p_hat}      # spatial dependence of f_000 at p_hat

# Derived quantities
v = p_hat/numpy.sqrt(p_hat**2 + 1)   # magnitude of the particle velocity
# Normalization of the distribution functionx
N = - numpy.sqrt(4 * numpy.pi/3) * 3 * u_one/v *  3*Q/(u_one * p_inj) \
    * r/(r - 1) * (p_hat/p_inj)**(-3*r/(r-1))

x_coords = inputs[0].Points[:,0]
outputArray = numpy.zeros(x_coords.shape[0])

# Upstream
outputArray[x_coords < 0] = N * numpy.exp(3 * u_one * nu_0 * x_coords[x_coords < 0]/(v*v))
# Downstream
outputArray[x_coords > 0] = 0

return outputArray""".format(
    Q=Q, p_inj=p_inj, r=r, u_one=u_one, nu_0=nu_0, p_hat=p_hat
)

python_calc_f_100_ana_x.MultilineExpression = python_calc_expression_f_100_x

# -----------------------
# Create a 'Plot Over Line'
# -----------------------

plotOverLine_f_x = PlotOverLine(
    registrationName="f(x) Plot", Input=python_calc_f_100_ana_x
)

# Extract min_x and max_x from the bounds
min_x = bounds[0]
max_x = bounds[1]

from math import log

ln_p_hat = log(p_hat)

# Specify where to plot over line
plotOverLine_f_x.Point1 = [min_x, ln_p_hat, 0.0]
plotOverLine_f_x.Point2 = [max_x, ln_p_hat, 0.0]

# Set SamplingPattern
plotOverLine_f_x.SamplingPattern = "Sample At Segment Centers"
# # The automatic ComputeTolerance does not work for large shock-grids,
# # so we set it manually
# plotOverLine_f_x.ComputeTolerance = 0
# # A tolerance of 1e-8 works well,
# # this is approx (dx_min / dx_max) / 1000,
# # with dx_min/max the minimum/maximum cell size in x
# plotOverLine_f_x.Tolerance = 1e-8

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
lineChartView_f_x.ViewSize = [width, height]
lineChartView_f_x.ViewTime = solution.TimestepValues[-1]  # last time step
# lineChartView_f_x.ChartTitle = "Spatial dependence"
lineChartView_f_x.LegendLocation = "TopLeft"
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
# SetActiveView(lineChartView_f_x)

# set active source
# SetActiveSource(plotOverLine_f_x)

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
fxPlotDisplay.SeriesLabel = [
    "f_000",
    "$f_{000}$",
    "f_100",
    "$f_{100}$",
    "f_000_ana_x",
    "$f_{000}$-ana",
    "f_100_ana_x",
    "$f_{100}$-ana",
]

# Adapt line thickness
fxPlotDisplay.SeriesLineThickness = [
    "f_000",
    "2",
    "f_000_ana_x",
    "2",
    "f_100",
    "2",
    "f_100_ana_x",
    "2",
]

# Line style
fxPlotDisplay.SeriesLineStyle = [
    "f_000",
    "1",
    "f_000_ana_x",
    "2",
    "f_100",
    "1",
    "f_100_ana_x",
    "2",
]  # ana dashed

# Color the plots
fxPlotDisplay.SeriesColor = [
    "f_000",
    "0.10980392156862745",
    "0.5843137254901961",
    "0.8039215686274",
    "f_000_ana_x",
    "0.25882352941176473",
    "0.23921568627450981",
    "0.6627450980392157",
    "f_100",
    "0.450980392157",
    "0.603921568627",
    "0.835294117647",
    "f_100_ana_x",
    "0.25882352941176473",
    "0.23921568627450981",
    "0.6627450980392157",
]

# Ensure that f_000 and f_000_ana_x are displayed
fxPlotDisplay.SeriesVisibility = [
    "f_000",
    "f_000_ana_x",
    "f_100",
    "f_100_ana_x",
]

# ----------------
# Save screenshot
# ----------------

print(f"Save screenshot '{results_folder}/spatial-distribution.png'")

# save screenshot
SaveScreenshot(
    filename=results_folder + "/spatial-distribution.png",
    location=vtkPVSession.DATA_SERVER,
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
    location=vtkPVSession.DATA_SERVER,
    ChooseArraysToWrite=2,
    PointDataArrays=["f_000", "f_000_ana_x"],
    Precision=6,
    UseScientificNotation=1,
)

# =============================================================================
# Plot over time (Time-dependent acceleration)
# =============================================================================

probeLocation = ProbeLocation(
    registrationName="PointEvaluation",
    Input=solution,
    ProbeType="Fixed Radius Point Source",
)

probeLocation.ProbeType.Center = [0.001, ln_p_hat, 0.0]

plotOverLine_f_t = PlotDataOverTime(
    registrationName="TemporalEvolution",
    OnlyReportSelectionStatistics=0,
    Input=probeLocation,
)
# # # ------------------------------------
# # Create new layout and LineChartView
# # ------------------------------------

# create new layout object 'f(t) Plot'
layout_f_t = CreateLayout(name="f(t) Plot")

# # Enter preview mode
layout_f_t.PreviewMode = [width, height]

# # layout/tab size in pixels
layout_f_t.SetSize(width, height)

# Create a new 'Line Chart View'
lineChartView_f_t = CreateView("XYChartView")
# NOTE: ViewSize needs to be set if pvbatch(pvpython) is used to produce the
# images
lineChartView_f_t.ViewSize = [width, height]
# lineChartView_f_t.ChartTitle = "Temporal dependence"
lineChartView_f_t.LegendLocation = "TopLeft"
lineChartView_f_t.LeftAxisTitle = "$f_{000}(t)$"
lineChartView_f_t.LeftAxisLogScale = 0
lineChartView_f_t.BottomAxisTitle = "$t$"
lineChartView_f_t.ChartTitleFontSize = 30
lineChartView_f_t.LeftAxisTitleFontSize = 24
lineChartView_f_t.BottomAxisTitleFontSize = 24
lineChartView_f_t.LegendFontSize = 18
lineChartView_f_t.LeftAxisLabelFontSize = 18
lineChartView_f_t.BottomAxisLabelFontSize = 18

# # assign view to a particular cell in the layout
AssignViewToLayout(view=lineChartView_f_t, layout=layout_f_t, hint=0)

# ---------------------
# Display PlotOverLine
# ---------------------

# show data in view
ftPlotDisplay = Show(
    plotOverLine_f_t, lineChartView_f_t, "XYChartRepresentation"
)

# Set which data to use for the x-axis
ftPlotDisplay.XArrayName = "Time"

# Set the labels of the plots
ftPlotDisplay.SeriesLabel = [
    "f_000 (id=0)",
    "$f_{000}$",
    "f_100 (id=0)",
    "$f_{100}$",
]

# Adapt line thickness
ftPlotDisplay.SeriesLineThickness = [
    "f_000 (id=0)",
    "2",
    "f_100 (id=0)",
]

ftPlotDisplay.SeriesColor = [
    "f_000 (id=0)",
    "0.10980392156862745",
    "0.5843137254901961",
    "0.8039215686274",
    "f_100 (id=0)",
    "0.450980392157",
    "0.603921568627",
    "0.835294117647",
]

ftPlotDisplay.SeriesVisibility = [
    "f_000 (id=0)",
    "f_100 (id=0)",
]

python_calc_f_000_ana_t = PythonCalculator(
    registrationName="f_000_ana_t", Input=plotOverLine_f_t
)
python_calc_f_000_ana_t.UseMultilineExpression = 1

python_calc_f_000_ana_t.MultilineExpression = """#########
# Physical Parameters
Q = {Q}              # number density of the injected particles
p_inj = {p_inj}      # magnitude of the injection momentum
r = {r}              # compression ratio of the shock
u_one = {u_one}      # shock velocity
nu_0 = {nu_0}        # inverse of the hall parameter h, i.e. h = omega_g/nu
p_hat = {p_hat}      # spatial dependence of f_000 at p_hat
delta_t = {delta_t}
# Normalization of the distribution functionx
N = numpy.sqrt(4 * numpy.pi) * 3*Q/(u_one * p_inj) \
    * r/(r - 1) * (p_hat/p_inj)**(-3*r/(r-1))

# Work around to get a numpy array with the time steps
# No idea if there is a better way
t_coords = numpy.zeros((1,Time.GetSize()), dtype=np.float64)
t_coords[:] = Time.GetArrays()
t_coords = t_coords.flatten()[1:] * delta_t # [1:] removes t = 0

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

return outputArray""".format(
    Q=Q, p_inj=p_inj, r=r, u_one=u_one, nu_0=nu_0, p_hat=p_hat, delta_t=delta_t
)
python_calc_f_000_ana_t.ArrayName = "f_000_ana_t"
PlotCalc = Show(python_calc_f_000_ana_t, lineChartView_f_t)
PlotCalc.SeriesLabel = ["f_000_ana_t (id=0)", "$f_{000}$-ana"]
PlotCalc.SeriesLineThickness = ["f_000_ana_t (id=0)", "2"]
PlotCalc.SeriesLineStyle = ["f_000_ana_t (id=0)", "2"]
PlotCalc.SeriesColor = [
    "f_000_ana_t (id=0)",
    "0.25882352941176473",
    "0.23921568627450981",
    "0.6627450980392157",
]

PlotCalc.SeriesVisibility = ["f_000_ana_t (id=0)"]

# ----------------
# Save screenshot
# ----------------

print(f"Save screenshot '{results_folder}/temporal-evolution.png'")

# save screenshot
SaveScreenshot(
    filename=results_folder + "/temporal-evolution.png",
    location=vtkPVSession.DATA_SERVER,
    viewOrLayout=layout_f_t,
)

# # ----------
# # Save data
# # ----------

print(f"Save data '{results_folder}/temporal-evolution.csv'")

# save data
SaveData(
    filename=results_folder + "/temporal-evolution.csv",
    proxy=plotOverLine_f_t,
    location=vtkPVSession.DATA_SERVER,
    Precision=6,
    UseScientificNotation=1,
)

# Render()
# Interact()
