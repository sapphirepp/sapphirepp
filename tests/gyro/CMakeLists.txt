add_test(NAME test_gyro_radius COMMAND ${CMAKE_BINARY_DIR}/examples/gyro_radius/gyro_radius ${CMAKE_BINARY_DIR}/examples/gyro_radius/parameter.prm gyro 0.1)
set_property(TEST test_gyro_radius PROPERTY LABELS gyro_radius fast_test)

add_test(NAME test_gyro_radius_mpi COMMAND mpirun ${CMAKE_BINARY_DIR}/examples/gyro_radius/gyro_radius ${CMAKE_BINARY_DIR}/examples/gyro_radius/parameter.prm gyro 0.1)
set_property(TEST test_gyro_radius_mpi PROPERTY LABELS gyro_radius fast_test mpi)

file(COPY mpi_performance.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY mpi_performance.csv DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME test_mpi_performance COMMAND bash mpi_performance.sh)
set_property(TEST test_mpi_performance PROPERTY LABELS gyro_radius slow_test mpi performance)

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
