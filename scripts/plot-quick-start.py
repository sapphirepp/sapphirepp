# trace generated using paraview version 5.13.2
# import paraview
# paraview.compatibility.major = 5
# paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
from paraview.util import *

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


# =============================================================================
# Plot 2D render view
# =============================================================================

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create new layout object '2D Plot'
# Using a new layout for every plot simplifies the scaling of fonts
layout1 = CreateLayout(name="2D Plot")

# Create a new 'Render View'
renderView1 = CreateView("RenderView")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=renderView1, layout=layout1, hint=0)

# set active view
SetActiveView(renderView1)

# show data in view
solutionDisplay = Show(solution, renderView1, "UnstructuredGridRepresentation")

# Enter preview mode
layout1.PreviewMode = [1024, 1024]

# layout/tab size in pixels
layout1.SetSize(1024, 1024)

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


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Exit preview mode
# layout1.PreviewMode = [0, 0]
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
