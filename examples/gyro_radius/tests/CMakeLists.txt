add_test(NAME test_gyro_radius COMMAND ../gyro_radius ../parameter.prm gyro 0.1)
set_property(TEST test_gyro_radius PROPERTY LABELS gyro_radius fast_test)

add_test(NAME test_gyro_radius_mpi COMMAND mpirun ../gyro_radius ../parameter.prm gyro 0.1)
set_property(TEST test_gyro_radius_mpi PROPERTY LABELS gyro_radius fast_test mpi)
