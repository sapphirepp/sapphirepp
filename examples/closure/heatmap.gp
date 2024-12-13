############################################################
# Plot a heatmap of the distribution function              #
############################################################
set term png enhanced size 1024,1024 font "STIX Two Math, 18"

# Define default values for input and output files
if (!exists("input_file")) input_file = 'results/closure/spherical_density_map_p_0_at_t_0.dat'
if (!exists("output_file")) output_file = 'results/closure/heatmap.png'

# Define output file from command-line argument or default value
set output output_file

set palette rgb 33,13,10 # rainbow
set view equal xyz
# set view 0,0
# Put coordinate system in the center of the sphere
set zeroaxis lt -1 dashtype 2 lc rgb '#E6000000'
set xyplane at 0
unset border
set xrange[-1.2:1.2]
set yrange[-1.2:1.2]
set zrange[-1.2:1.2]

# set xtics axis nomirror front
set label 1 'V_x' at first 1.3,0,0 front center
set label 2 'V_y' at first 0,0,1.3 front center
set label 3 'V_z' at first 0,-1.22,0 front center

unset xtics
unset ytics
unset ztics

unset key # Removes the file names from the plot
set pm3d implicit at s
set pm3d depthorder
set linetype 10 lc rgb '#E6000000' # first two hex determines the
set pm3d interpolate 1,1 flush begin noftriangles border lt 10

set colorbox vertical user origin 0.,0.3 size 0.02,0.4
# Change range of colormap
set cbrange[-0.05:0.25]
set cbtics scale 0.5
set cbtics (-0.05, 0., 0.05, 0.1, 0.15, 0.2) out nomirror

# Define input file from command-line argument or default value
splot input_file with lines

unset colorbox