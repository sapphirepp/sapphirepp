# for mpi_processes in 1 2 4 6 8 #TODO: Something weird happends for -np 4, 6
for mpi_processes in 1 2 8 12
do
    mpirun -np $mpi_processes ../gyro_radius ../parameter.prm mpi_processes mpi_performance.csv
done