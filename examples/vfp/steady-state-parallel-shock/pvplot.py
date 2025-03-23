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
# Physical Parameters for the analytical solution
# =============================================================================
Q = 1.              # number density of the injected particles
p_inj = 2.           # magnitude of the injection momentum
r = 4                # compression ratio of the shock
u_one = 0.0167       # shock velocity
nu_0 = 1             # inverse of the hall parameter h, i.e. h = omega_g/nu
p_hat = 50           # spatial dependence of f_000 is shown at p = p_hat

# =============================================================================
# Properties of the images of the plots
# =============================================================================
width = 1280
height = 720

# =============================================================================
# Load pvtu files
# =============================================================================

from paraview.util import *

pvtu_files = paraview.util.Glob(results_folder + "/f_*.pvtu")
if not pvtu_files:
    raise FileNotFoundError(
        f"No .pvtu files found matching '{results_folder}/f_*.pvtu'"
    )

# create a new 'XML Partitioned Unstructured Grid Reader'
solution = XMLPartitionedUnstructuredGridReader(
    registrationName="f",
    FileName=pvtu_files,
)

# for pvtu, the time variable does not work
solution.TimeArray = "None"
# Only access 'f_000'
solution.PointArrayStatus = ["f_000"]

# Update Pipeline
solution.UpdatePipeline()

# Get the bounds
# Get data information from the solution
solution_data = servermanager.Fetch(solution)
bounds = solution_data.GetBounds()

# =============================================================================
# Plot f(x)
# =============================================================================

# -----------------------
# Compute the analytical solution
# -----------------------

python_calc_ana_x = PythonCalculator(registrationName="f_000_ana_x", Input=solution)
python_calc_ana_x.ArrayName = 'f_000_ana_x'
python_calc_ana_x.UseMultilineExpression = 1

python_cal_expression_x = """#########
# Physical Parameters
Q = {Q}              # number density of the injected particles
p_inj = {p_inj}      # magnitude of the injection momentum
r = {r}              # compression ratio of the shock
u_one = {u_one}      # shock velocity
nu_0 = {nu_0}        # inverse of the hall parameter h, i.e. h = omega_g/nu
p_hat = {p_hat}      # spatial dependence of f_000 at p_hat

# Derived quantities
v = p_hat/numpy.sqrt(p_hat**2 + 1)   # magnitude of the particle velocity
nu = nu_0/p_hat                      # scattering frequency for p_hat
# Normalization of the distribution function
N = 3*Q/(numpy.sqrt(4 * numpy.pi) * u_one * p_inj**3) \
    * r/(r - 1) * (p_hat/p_inj)**(-3*r/(r-1))  

x_coords = inputs[0].Points[:,0]
outputArray = numpy.zeros(x_coords.shape[0])

# Upstream
outputArray[x_coords < 0] = N * numpy.exp(3 * u_one * nu * x_coords[x_coords < 0]/(v*v))
# Downstream
outputArray[x_coords > 0] = N

return outputArray""".format(Q = Q,
                             p_inj = p_inj,
                             r = r,
                             u_one = u_one,
                             nu_0 = nu_0,
                             p_hat = p_hat)

python_calc_ana_x.MultilineExpression = python_cal_expression_x

# -----------------------
# Create a 'Plot Over Line'
# -----------------------

plotOverLine_f_x = PlotOverLine(registrationName="f(x) Plot", Input=python_calc_ana_x)

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


# =============================================================================
# Plot f(p)
# =============================================================================

# -----------------------
# Compute the analytical solution
# -----------------------

python_calc_ana_p = PythonCalculator(registrationName="f_000_ana_p", Input=solution)
python_calc_ana_p.ArrayName = 'f_000_ana_p'
python_calc_ana_p.UseMultilineExpression = 1

python_cal_expression_p = """#########
# Physical Parameters
Q = {Q}              # number density of the injected particles
p_inj = {p_inj}      # magnitude of the injection momentum
r = {r}              # compression ratio of the shock
u_one = {u_one}      # shock velocity

y_coords = numpy.exp(inputs[0].Points[:,1])
outputArray = 3*Q/(numpy.sqrt(4 * numpy.pi) * u_one * p_inj**3) \
    * r/(r - 1) * (y_coords/p_inj)**(-3*r/(r-1)) 

return outputArray""".format(Q = Q,
                             p_inj = p_inj,
                             r = r,
                             u_one = u_one)

python_calc_ana_p.MultilineExpression = python_cal_expression_p

# -----------------------
# Create a 'Plot Over Line'
# -----------------------

# create a new 'Plot Over Line'
plotOverLine_f_p = PlotOverLine(
    registrationName="f(p) Plot", Input=python_calc_ana_p)

# Get the bounds in y
min_y = bounds[2]
max_y = bounds[3]

# Specify where to plot over line
# Properties modified on plotOverLine_f_p
plotOverLine_f_p.Point1 = [0.1, min_y, 0.0]
plotOverLine_f_p.Point2 = [0.1, max_y, 0.0]

# Set SamplingPattern
plotOverLine_f_p.SamplingPattern = "Sample At Segment Centers"

# ------------------------------------
# Create new layout and LineChartView
# ------------------------------------

# create new layout object 'f(p) Plot'
layout_f_p = CreateLayout(name="f(p) Plot")

# Enter preview mode
layout_f_p.PreviewMode = [width, height]

# layout/tab size in pixels
layout_f_p.SetSize(width, height)

# Create a new 'Line Chart View'
lineChartView_f_p = CreateView("XYChartView")
# NOTE: ViewSize needs to be set if pvbatch(pvpython) is used to produce the
# images
lineChartView_f_p.ViewSize = [width,height] 
# lineChartView_f_p.ChartTitle = "f(p) Plot"
lineChartView_f_p.LeftAxisTitle = "$f_{000}(\\ln p)$"
lineChartView_f_p.BottomAxisTitle = "$\\ln p$"
lineChartView_f_p.ChartTitleFontSize = 30
lineChartView_f_p.LeftAxisTitleFontSize = 24
lineChartView_f_p.BottomAxisTitleFontSize = 24
lineChartView_f_p.LegendFontSize = 18
lineChartView_f_p.LeftAxisLabelFontSize = 18
lineChartView_f_p.BottomAxisLabelFontSize = 18

lineChartView_f_p.LeftAxisUseCustomRange = 1
lineChartView_f_p.LeftAxisRangeMinimum = 1e-6
lineChartView_f_p.LeftAxisRangeMaximum = 14.5
lineChartView_f_p.LeftAxisLogScale = 1 # Use logarithmic y axis

lineChartView_f_p.BottomAxisUseCustomRange = 1
lineChartView_f_p.BottomAxisRangeMinimum = 0.
lineChartView_f_p.BottomAxisRangeMaximum = 4.8

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

# Properties modified on fpPlotDisplay
fpPlotDisplay.XArrayName = "Points_Y"
# Set the labels of the plots
fpPlotDisplay.SeriesLabel =  [ 'f_000', '$f_{000}$', 'f_000_ana_p', '$f_{000}$-ana']

# Adapt line thickness
fpPlotDisplay.SeriesLineThickness = ['f_000', '2', 'f_000_ana_p', '2']

# Line style
fpPlotDisplay.SeriesLineStyle = ['f_000', '1', 'f_000_ana_p', '2'] # ana dashed

# Color the plots
fpPlotDisplay.SeriesColor = [
    'f_000', '0.10980392156862745', '0.5843137254901961', '0.8039215686274',
    'f_000_ana_p', '0.25882352941176473', '0.23921568627450981', '0.6627450980392157'
]

# Ensure that f_000 and f_000_ana_p are displayed
fpPlotDisplay.SeriesVisibility = ['f_000', 'f_000_ana_p']

# ----------------
# Save screenshot
# ----------------

print(f"Save screenshot '{results_folder}/particle-spectrum.png'")

# save screenshot
SaveScreenshot(
    filename=results_folder + "/particle-spectrum.png",
    viewOrLayout=layout_f_p,
    location = vtkPVSession.DATA_SERVER,
)


# ----------
# Save data
# ----------

print(f"Save data '{results_folder}/particle-spectrum.csv'")

# save data
SaveData(
    filename=results_folder + "/particle-spectrum.csv",
    proxy=plotOverLine_f_p,
    location = vtkPVSession.DATA_SERVER,
    ChooseArraysToWrite=2,
    PointDataArrays=['f_000', 'f_000_ana_p'],
    Precision=6,
    UseScientificNotation=1,
)

# =============================================================================
# Plot p^4 f(p) if a subfolder 'scaled' exits in the results_folder
# =============================================================================
from os.path import isdir
if isdir(results_folder + '/scaled'):
    # ==========================================================================
    # Load pvtu files
    # ==========================================================================

    pvtu_files_scaled = paraview.util.Glob(results_folder + "/scaled/g_*.pvtu")
    if not pvtu_files:
        raise FileNotFoundError(
            f"No .pvtu files found matching '{results_folder}/scaled/g_*.pvtu'"
        )

    # create a new 'XML Partitioned Unstructured Grid Reader'
    solution_scaled = XMLPartitionedUnstructuredGridReader(
        registrationName="g_scaled",
        FileName=pvtu_files_scaled,
    )

    # for pvtu, the time variable does not work
    solution_scaled.TimeArray = "None"
    # Only access 'f_000'
    solution_scaled.PointArrayStatus = ["g_000"]

    # Get the bounds
    # Get data information from the solution
    solution_scaled_data = servermanager.Fetch(solution_scaled)
    bounds_scaled = solution_scaled_data.GetBounds()

    # -----------------------
    # Scale f_000, i.e. multiply with p to get p^4 f_000
    # -----------------------
    calc_scale_f_000 = Calculator(
        registrationName='g_000_scaled',
        Input=solution_scaled
    )
    calc_scale_f_000.ResultArrayName = 'g_000_scaled'
    calc_scale_f_000.Function = 'exp(coordsY) * g_000'
    
    # -----------------------
    # Compute the analytical solution
    # -----------------------
    python_calc_ana_p_scaled = PythonCalculator(
        registrationName="g_000_ana_p_scaled",
        Input=calc_scale_f_000
    )
    python_calc_ana_p_scaled.ArrayName = 'g_000_ana_p_scaled'
    python_calc_ana_p_scaled.UseMultilineExpression = 1

    python_cal_expression_p_scaled = """#########
# Physical Parameters
Q = {Q}               # number density of the injected particles
p_inj = {p_inj}       # magnitude of the injection momentum
r = {r}               # compression ratio of the shock
u_one = {u_one}       # shock velocity

N = 3*Q/(numpy.sqrt(4 * numpy.pi) * u_one * p_inj**3) \
    * r/(r - 1) * (p_inj)**(3*r/(r-1)) 

outputArray = N * numpy.ones(inputs[0].PointData['g_000'].shape[0])

return outputArray""".format(Q = Q,
                             p_inj = p_inj,
                             r = r,
                             u_one = u_one)
    
    python_calc_ana_p_scaled.MultilineExpression = python_cal_expression_p_scaled

    # -----------------------
    # Create a 'Plot Over Line'
    # -----------------------

    # create a new 'Plot Over Line'
    plotOverLine_f_p_scaled = PlotOverLine(
        registrationName="f(p) Plot (Scaled)", Input=python_calc_ana_p_scaled)

    # Get the bounds in y
    min_y = bounds_scaled[2]
    max_y = bounds_scaled[3]

    # Specify where to plot over line
    # Properties modified on plotOverLine_f_p
    plotOverLine_f_p_scaled.Point1 = [0.1, min_y, 0.0]
    plotOverLine_f_p_scaled.Point2 = [0.1, max_y, 0.0]

    # Set SamplingPattern
    plotOverLine_f_p_scaled.SamplingPattern = "Sample At Segment Centers"


    # ------------------------------------
    # Create new layout and LineChartView
    # ------------------------------------

    # create new layout object 'f(p) Plot'
    layout_f_p_scaled = CreateLayout(name="f(p) Plot (Scaled)")

    # Enter preview mode
    layout_f_p_scaled.PreviewMode = [width, height]

    # layout/tab size in pixels
    layout_f_p_scaled.SetSize(width, height)

    # Create a new 'Line Chart View'
    lineChartView_f_p_scaled = CreateView("XYChartView")
    # NOTE: ViewSize needs to be set if pvbatch(pvpython) is used to produce the
    # images
    lineChartView_f_p_scaled.ViewSize = [width, height] 
    # lineChartView_f_p_scaled.ChartTitle = "f(p) Plot"
    lineChartView_f_p_scaled.BottomAxisTitle = "$\\ln p$"
    lineChartView_f_p_scaled.ChartTitleFontSize = 30
    lineChartView_f_p_scaled.LeftAxisTitleFontSize = 24
    lineChartView_f_p_scaled.BottomAxisTitleFontSize = 24
    lineChartView_f_p_scaled.LegendLocation = 'BottomLeft'
    lineChartView_f_p_scaled.LegendFontSize = 18
    lineChartView_f_p_scaled.LeftAxisLabelFontSize = 18
    lineChartView_f_p_scaled.BottomAxisLabelFontSize = 18

    lineChartView_f_p_scaled.LeftAxisUseCustomRange = 1
    lineChartView_f_p_scaled.LeftAxisRangeMinimum = 0.01
    lineChartView_f_p_scaled.LeftAxisRangeMaximum = 140
    lineChartView_f_p_scaled.LeftAxisLogScale = 0 # Use logarithmic y axis
    
    lineChartView_f_p_scaled.BottomAxisUseCustomRange = 1
    lineChartView_f_p_scaled.BottomAxisRangeMinimum = 0.
    lineChartView_f_p_scaled.BottomAxisRangeMaximum = 4.8

    # assign view to a particular cell in the layout
    AssignViewToLayout(view=lineChartView_f_p_scaled,
                       layout=layout_f_p_scaled,
                       hint=0)
    
    # set active view
    SetActiveView(lineChartView_f_p_scaled)
    
    # set active source
    SetActiveSource(plotOverLine_f_p_scaled)


    # ---------------------
    # Display PlotOverLine
    # ---------------------
    
    # show data in view
    fpScaledPlotDisplay = Show(
        plotOverLine_f_p_scaled, lineChartView_f_p_scaled, "XYChartRepresentation"
    )
    
    # Properties modified on fpScaledPlotDisplay
    fpScaledPlotDisplay.XArrayName = "Points_Y"
    # Set the labels of the plots
    fpScaledPlotDisplay.SeriesLabel =  [ 'g_000', '$p^{3}f_{000}$',
                                         'g_000_scaled','$p^{4}f_{000}$',
                                         'g_000_ana_p_scaled', '$p^{4}f_{000}$-ana']
    
    # Adapt line thickness
    fpScaledPlotDisplay.SeriesLineThickness = ['g_000', '2',
                                               'g_000_scaled', '2',
                                               'g_000_ana_p_scaled', '2']
    
    # Line style
    fpScaledPlotDisplay.SeriesLineStyle = ['g_000', '1',
                                           'g_000_scaled', '1',
                                           'g_000_ana_p_scaled', '2'] # ana dashed
    
    # Color the plots
    fpScaledPlotDisplay.SeriesColor = [
        'g_000', '0.10980392156862745', '0.5843137254901961', '0.8039215686274',
        'g_000_scaled', '0.3058823529411765', '0.8509803921568627', '0.9176470588235294',
        'g_000_ana_p_scaled', '0.25882352941176473', '0.23921568627450981', '0.6627450980392157'
    ]
    
    # Ensure that f_000 and f_000_ana_p are displayed
    fpScaledPlotDisplay.SeriesVisibility = ['g_000', 'g_000_scaled', 'g_000_ana_p_scaled']
    
    # ----------------
    # Save screenshot
    # ----------------
    
    print(f"Save screenshot '{results_folder}/scaled-particle-spectrum.png'")
    
    # save screenshot
    SaveScreenshot(
        filename=results_folder + "/scaled-particle-spectrum.png",
        viewOrLayout=layout_f_p_scaled,
        location = vtkPVSession.DATA_SERVER,
    )
    
    
    # ----------
    # Save data
    # ----------
    
    print(f"Save data '{results_folder}/scaled-particle-spectrum.csv'")
    
    # save data
    SaveData(
        filename=results_folder + "/scaled-particle-spectrum.csv",
        proxy=plotOverLine_f_p_scaled,
        location = vtkPVSession.DATA_SERVER,
        ChooseArraysToWrite=3,
        PointDataArrays=['g_000', 'g_000_scaled', 'g_000_ana_p_scaled'],
        Precision=6,
        UseScientificNotation=1,
    )

## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------
