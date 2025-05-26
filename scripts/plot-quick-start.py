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

# Update Pipeline
solution.UpdatePipelineInformation()


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
layout2D.PreviewMode = [1024, 1024]
# layout/tab size in pixels
layout2D.SetSize(1024, 1024)

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


# ------------------
# Show DataAxesGrid
# ------------------

solutionDisplay.DataAxesGrid.GridAxesVisibility = 1
solutionDisplay.DataAxesGrid.XTitle = "$x$"
solutionDisplay.DataAxesGrid.YTitle = "$\\ln p$"
# Only show Axes Min-X//Y/Z
solutionDisplay.DataAxesGrid.AxesToLabel = 7
# solutionDisplay.DataAxesGrid.FacesToRender = 7
# Set default font size: 24 for title, 18 for label
solutionDisplay.DataAxesGrid.XTitleFontSize = 24
solutionDisplay.DataAxesGrid.YTitleFontSize = 24
solutionDisplay.DataAxesGrid.XLabelFontSize = 18
solutionDisplay.DataAxesGrid.YLabelFontSize = 18
# Use gray color for label for good visibility in both light and dark mode
solutionDisplay.DataAxesGrid.XTitleColor = [0.5, 0.5, 0.5]
solutionDisplay.DataAxesGrid.YTitleColor = [0.5, 0.5, 0.5]
solutionDisplay.DataAxesGrid.XLabelColor = [0.5, 0.5, 0.5]
solutionDisplay.DataAxesGrid.YLabelColor = [0.5, 0.5, 0.5]
solutionDisplay.DataAxesGrid.GridColor = [0.5, 0.5, 0.5]


# update the view to ensure updated data information
renderView1.Update()


# ----------------------
# Create color bar plot
# ----------------------

# set scalar coloring
ColorBy(solutionDisplay, ("POINTS", "f_000"))

# show color bar/color legend
solutionDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'f_000'
f_000LUT = GetColorTransferFunction("f_000")
# get opacity transfer function/opacity map for 'f_000'
f_000PWF = GetOpacityTransferFunction("f_000")
# get 2D transfer function for 'f_000'
f_000TF2D = GetTransferFunction2D("f_000")

# Rescale transfer function
f_000LUT.RescaleTransferFunction(1e-2, 10.0)
# Rescale transfer function
f_000PWF.RescaleTransferFunction(1e-2, 10.0)
# Rescale 2D transfer function
f_000TF2D.RescaleTransferFunction(1e-2, 10.0, 0.0, 1.0)

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
f_000LUTColorBar.Position = [0.15, 0.55]


# Go to last simulation time step
animationScene1.GoToLast()


# ----------------
# Save screenshot
# ----------------

# Save screenshot in results folder (at data server).
print(f"Save screenshot '{results_folder}/quick-start-2D.png'")

# save screenshot
SaveScreenshot(
    filename=results_folder + "/quick-start-2D.png",
    viewOrLayout=renderView1,
    location=vtkPVSession.DATA_SERVER,
    TransparentBackground=1,
)


# ----------------
# Save Animation
# ----------------

print(f"Save animation '{results_folder}/quick-start-2D.*.png'")

# save animation
SaveAnimation(
    filename=results_folder + "/quick-start-2D.png",
    viewOrLayout=renderView1,
    location=vtkPVSession.DATA_SERVER,
    TransparentBackground=1,
    FrameStride=1,
)


# =============================================================================
# Plot f(x)
# =============================================================================

# create a new 'Plot Over Line'
plotOverLine_f_x = PlotOverLine(registrationName="f(x) Plot", Input=solution)

# Get the bounds in x
# Fetch data information from the solution
solution_data = servermanager.Fetch(solution)
# Get bounds of the data
bounds = solution_data.GetBounds()
# Extract min_x and max_x from the bounds
min_x = bounds[0]
max_x = bounds[1]

# Properties modified on plotOverLine_f_x
plotOverLine_f_x.Point1 = [min_x, 0.05, 0.0]
plotOverLine_f_x.Point2 = [max_x, 0.05, 0.0]


# -----------------------
# Change SamplingPattern
# -----------------------

# Default
plotOverLine_f_x.SamplingPattern = "Sample Uniformly"
plotOverLine_f_x.Resolution = 1000

# Sample at cell boundaries for full DG data (preferred)
plotOverLine_f_x.SamplingPattern = "Sample At Cell Boundaries"

# Sample at cell centers for smooth line
plotOverLine_f_x.SamplingPattern = "Sample At Segment Centers"


# ------------------------------------
# Create new layout and LineChartView
# ------------------------------------

# create new layout object 'f(x) Plot'
layout_f_x = CreateLayout(name="f(x) Plot")

# Create a new 'Line Chart View'
lineChartView_f_x = CreateView("XYChartView")
# lineChartView_f_x.ChartTitle = "$f(x)$ Plot"
lineChartView_f_x.LeftAxisTitle = "$f_{lms}(x)$"
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

# Enter preview mode
layout_f_x.PreviewMode = [1280, 720]
# layout/tab size in pixels
layout_f_x.SetSize(1280, 720)

# Properties modified on fxPlotDisplay
fxPlotDisplay.XArrayName = "Points_X"
fxPlotDisplay.SeriesLineThickness = [
    "f_000",
    "3",
    "f_100",
    "3",
]


# ----------------
# Save screenshot
# ----------------

print(f"Save screenshot '{results_folder}/quick-start-f-x.png'")

# save screenshot
SaveScreenshot(
    filename=results_folder + "/quick-start-f-x.png",
    viewOrLayout=layout_f_x,
    location=vtkPVSession.DATA_SERVER,
    TransparentBackground=1,
)


# ----------
# Save data
# ----------

print(f"Save data '{results_folder}/quick-start-f-x.csv'")

# save data
SaveData(
    filename=results_folder + "/quick-start-f-x.csv",
    proxy=plotOverLine_f_x,
    location=vtkPVSession.DATA_SERVER,
    ChooseArraysToWrite=1,
    PointDataArrays=["f_000", "f_100"],
    Precision=5,
    UseScientificNotation=1,
)


# =============================================================================
# Create Calculator p^4 f_000
# =============================================================================

# create a new 'Calculator'
calculator_p4f = Calculator(registrationName="p^4 f_000", Input=solution)

# Properties modified on calculator_p4f
calculator_p4f.ResultArrayName = "$p^4 f_{000}$"
calculator_p4f.Function = "exp(4*coordsY) * f_000"


# =============================================================================
# Plot p^4 f(p)
# =============================================================================

# create a new 'Plot Over Line'
plotOverLine_f_p = PlotOverLine(
    registrationName="f(p) Plot", Input=calculator_p4f
)

# Get the bounds in y
# Fetch data information from the solution
solution_data = servermanager.Fetch(solution)
# Get bounds of the data
bounds = solution_data.GetBounds()
# Extract min_x and max_x from the bounds
min_y = bounds[2]
max_y = bounds[3]

# Properties modified on plotOverLine_f_p
plotOverLine_f_p.Point1 = [0.1, min_y, 0.0]
plotOverLine_f_p.Point2 = [0.1, max_y, 0.0]
plotOverLine_f_p.SamplingPattern = "Sample At Segment Centers"


# ------------------------------------
# Create new layout and LineChartView
# ------------------------------------

# create new layout object 'f(p) Plot'
layout_f_p = CreateLayout(name="f(p) Plot")

# Create a new 'Line Chart View'
lineChartView_f_p = CreateView("XYChartView")
# lineChartView_f_p.ChartTitle = "f(p) Plot"
lineChartView_f_p.LeftAxisTitle = "$p^4 f_{000}(\\ln p)$"
lineChartView_f_p.BottomAxisTitle = "$\\ln p$"
lineChartView_f_p.ChartTitleFontSize = 30
lineChartView_f_p.LeftAxisTitleFontSize = 24
lineChartView_f_p.BottomAxisTitleFontSize = 24
lineChartView_f_p.LegendFontSize = 18
lineChartView_f_p.LeftAxisLabelFontSize = 18
lineChartView_f_p.BottomAxisLabelFontSize = 18

# assign view to a particular cell in the layout
AssignViewToLayout(view=lineChartView_f_p, layout=layout_f_p, hint=0)

# set active view
SetActiveView(lineChartView_f_p)

# set active source
SetActiveSource(plotOverLine_f_p)


# ---------------------
# Display PlotOverLine
# ---------------------

# show data in view
fpPlotDisplay = Show(
    plotOverLine_f_p, lineChartView_f_p, "XYChartRepresentation"
)

# Enter preview mode
layout_f_p.PreviewMode = [1280, 720]
# layout/tab size in pixels
layout_f_p.SetSize(1280, 720)

# Properties modified on fpPlotDisplay
fpPlotDisplay.XArrayName = "Points_Y"
fpPlotDisplay.SeriesVisibility = ["$p^4 f_{000}$"]
fpPlotDisplay.SeriesLineThickness = [
    "$p^4 f_{000}$",
    "3",
]


# -----------------------
# Use logarithmic y axis
# -----------------------

# Properties modified on lineChartView_f_p
lineChartView_f_p.LeftAxisUseCustomRange = 1
lineChartView_f_p.LeftAxisRangeMaximum = 16.0
lineChartView_f_p.LeftAxisRangeMinimum = 1e-2
lineChartView_f_p.LeftAxisLogScale = 1


# ----------------
# Save screenshot
# ----------------

print(f"Save screenshot '{results_folder}/quick-start-f-p.png'")

# save screenshot
SaveScreenshot(
    filename=results_folder + "/quick-start-f-p.png",
    viewOrLayout=layout_f_p,
    location=vtkPVSession.DATA_SERVER,
    TransparentBackground=1,
)


# ----------
# Save data
# ----------

print(f"Save data '{results_folder}/quick-start-f-p.csv'")

# save data
SaveData(
    filename=results_folder + "/quick-start-f-p.csv",
    proxy=plotOverLine_f_p,
    location=vtkPVSession.DATA_SERVER,
    ChooseArraysToWrite=1,
    PointDataArrays=["f_000", "f_100", "$p^4 f_{000}$"],
    Precision=5,
    UseScientificNotation=1,
)


# =============================================================================
# ParaView Python and NumPy Introduction: Calculate Spectral Index
# =============================================================================

import numpy as np
from paraview.vtk.util import numpy_support

# Fetch the data from the plotOverLine object
plot_over_line_data = servermanager.Fetch(plotOverLine_f_p)

# Get the points and f_000 array
points_vtk = plot_over_line_data.GetPoints()
f_000_vtk = plot_over_line_data.GetPointData().GetArray("f_000")

# Convert points and f_000 to numpy array
points = numpy_support.vtk_to_numpy(points_vtk.GetData())
f_000 = numpy_support.vtk_to_numpy(f_000_vtk)

# Extract ln_p array array from y-coordinates
ln_p = points[:, 1]

# Filter out data below 2x injection momentum (p_inj = 1)
mask = ln_p > np.log(2.0)
ln_p = ln_p[mask]
f_000 = f_000[mask]

# Calculate ln(f)
ln_f = np.log(f_000)

# Calculate log-log-slope of spectrum to find the spectral index
spectral_index = (ln_f[-1] - ln_f[0]) / (ln_p[-1] - ln_p[0])

print(f"Spectral Index: s = {spectral_index}")


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Exit preview mode
# layout2D.PreviewMode = [0, 0]
# layout_f_x.PreviewMode = [0, 0]
# layout_f_p.PreviewMode = [0, 0]
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------
