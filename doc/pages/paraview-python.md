# ParaView Python Introduction {#paraview-python}

Often, we need to produce the same plots for multiple datasets
or perform additional manipulations on the data,
such as calculating the spectral index.
A common approach is to use Python for these tasks,
often in combination with `matplotlib` for visualization.
Unfortunately, the output of @sapphire is not straightforward to manipulate.
See the discussion on discontinuous Galerkin on the [next page](#visualization-tips).
Therefore, we discourage directly reading the @sapphire results in Python.
Instead, we recommend using tools like @paraview
and [VisIt](https://visit-dav.github.io/visit-website/)
that can handle these subtleties.
However, relying on a graphical user interface (GUI)
to create repetitive plots
or perform detailed analysis
can be cumbersome and disrupt existing workflows.
Fortunately, @paraview comes with a Python interface
that can bridge this gap.

The @paraview Python interface has two main aspects:

1. Executing Python scripts directly in the @paraview application,
   e.g., as macros to simplify repetitive tasks.
2. Running scripts from a shell
   using `pvpython` or `pvbatch`
   to create plots without a GUI.
   - `pvpython` is the Python interpreter that runs @paraview’s Python scripts.
     It links the ParaView Python package to the interpreter.
   - `pvbatch` is similar,
     but while `pvpython` is meant for interactive scripts,
     `pvbatch` is designed for batch processing.
     Additionally, when running on computing resources with MPI capabilities,
     `pvbatch` can be run in parallel.

Both executables come installed with @paraview.
Unfortunately, including other Python packages with these interpreters is not straightforward.
([See here](https://www.paraview.org/paraview-docs/latest/python/quick-start.html)
for how to use the @paraview Python package outside `pvpython`.)
Using an installation with
[conda/conda-forge](#paraview-installation)
simplifies this
by linking @paraview and the Python environment.

To get started with @paraview Python,
we provide some scripts as a starting point.
These scripts mainly serve to reproduce the plots in the [examples](#examples)
but also introduce some advanced concepts.
We illustrate their usage using the [quick-start](#quick-start) example.
The script for the quick-start is located at
[scripts/plot-quick-start.py](https://github.com/sapphirepp/sapphirepp/blob/main/scripts/plot-quick-start.py).
This script takes the path to the results,
i.e., `/path/to/sapphirepp/results/01`,
either as a command line argument
or as prompted input from the user.
It loads the data from that folder,
reproduces the plots from the
[ParaView Tutorial](#paraview-tutorial),
and saves them as `.png` files in the same folder.
This script can be executed in multiple ways:

- **Execute as Script in ParaView:**
  1. Open the *Python Shell* in @paraview
     via **View** > **Python Shell**.
  2. Click on **Run Script**
     and select `sapphirepp/scripts/plot-quick-start.py`.
  3. When prompted, enter the path to the results folder,
     e.g., `/path/to/sapphirepp/results/01`.
  4. The plots will open automatically in @paraview.

- **Open in Script Editor in ParaView:**
  1. Open the *Script Editor*
     via **Tools** > **Python Script Editor**.
  2. Open the script by selecting **File** > **Open**
     and choosing `sapphirepp/scripts/plot-quick-start.py`.
  3. Execute the script via **File** > **Run**.
  4. When prompted, enter the path to the results folder,
     and the plots will be created automatically.

- **Import as Macro in ParaView:**
  1. Go to **Macros** > **Import new macro...**
     and select the script.
  2. Execute the macro via **Macros** > **visualize-quick-start**.
  3. When prompted, enter the path to the results folder,
     and the plots will be created automatically.

- **Execute with `pvpython`:**
  1. Open a terminal shell
     and ensure `pvpython` is installed and in `PATH`,
     e.g., by using `conda activate ParaView`.
  2. Execute the script
     and provide the results path as a command line argument:

     ```shell
     pvpython scripts/plot-quick-start.py results/01
     ```

  3. This opens a temporary @paraview preview window,
     which closes after the script finishes.
  4. Note, executing `pvpython` for the first time can take a bit longer,
     as a @paraview instance has to be started in the background.

- **Execute with `pvbatch`:**
  1. Open a terminal shell
     and ensure `pvpython` is installed and in `PATH`,
     e.g., by using `conda activate ParaView`.
  2. Execute the script
     and provide the results path as a command line argument:

     ```shell
     pvbatch scripts/plot-quick-start.py results/01
     ```

  3. This does not open a @paraview preview window.

@note The first three options also work while connected to a
      [Remote ParaView Server](https://docs.paraview.org/en/latest/ReferenceManual/parallelDataVisualization.html).
  
The script creates a 2D plot `quick-start-2D.png`,
as well as an animated time series,
`quick-start-2D.0000-0200.png`.
To convert this into an animated `.gif`,
you can use
[ImageMagick Convert](https://imagemagick.org/script/convert.php):

```shell
magick -delay 1 -loop 0 results/01/quick-start-2D.*.png results/01/quick-start-2D.gif
rm results/01/quick-start-2D.*.png
```

The f(x) and f(p) plots created in the [ParaView tutorial](#paraview-tutorial)
are saved as `quick-start-f-x.png` and `quick-start-f-p.png`, respectively.
Additionally, the underlying `PlotOverLine` data is saved as `csv` files,
`quick-start-f-x.csv` and `quick-start-f-p.csv`.
This 1D data can be used for further analysis,
as @paraview has already interpolated the data
and tabulated it in an easy-to-understand format.
However, be aware of the quirks of the `PlotOverLine` feature,
as discussed on the [next page](#paraview-visualization).

Last, the script also calculates the spectral index $s$:
`Spectral Index: s = -4.049119472503662`.
This section in the
[plot-quick-start.py](https://github.com/sapphirepp/sapphirepp/blob/main/scripts/plot-quick-start.py)
script serves as an introduction
on how to combine [NumPy](https://numpy.org/) with @paraview.
It utilizes the
[`vtk.util.numpy_support`](https://docs.vtk.org/en/latest/api/python/vtkmodules/vtkmodules.util.numpy_support.html)
package
to facilitate the conversion between VTK and NumPy arrays.

To further familiarize yourself with @paraview Python,
we recommend using the **Tools** > **Start Trace** option,
which generates Python scripts from your interactions with the GUI.
Additionally, we refer to the relevant chapters in the
[ParaView User's Guide](https://docs.paraview.org/en/latest/UsersGuide/introduction.html#getting-started-with-pvpython)
and the
[ParaView Reference Manual](https://docs.paraview.org/en/latest/ReferenceManual/parallelDataVisualization.html#sec-usingpvbatch),
as well as the following tutorials:
[Batch Python Scripting](https://docs.paraview.org/en/latest/Tutorials/SelfDirectedTutorial/batchPythonScripting.html)
and
[Python & Batch: ParaView & Python](https://docs.paraview.org/en/latest/Tutorials/ClassroomTutorials/pythonAndBatchParaViewAndPython.html).
The documentation for the Python module can be found here:
[`paraview.simple` documentation](https://www.paraview.org/paraview-docs/latest/python/paraview.simple.html).

<div class="section_buttons">

| Previous                                |                                                     Next |
|:----------------------------------------|---------------------------------------------------------:|
| [ParaView Tutorial](#paraview-tutorial) | [Tips and Tricks for Visualization](#visualization-tips) |

</div>

---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2025-02-19
