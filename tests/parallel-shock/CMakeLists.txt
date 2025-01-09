file(COPY parameter.prm
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)


add_executable(test-parallel-shock test-parallel-shock.cpp ${VFP_SOURCES})
deal_ii_setup_target(test-parallel-shock)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/parallel-shock)
target_include_directories(test-parallel-shock PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp
)
target_include_directories(test-parallel-shock PUBLIC ${CONFIG_HEADER_DIR})
target_link_libraries(test-parallel-shock UtilsLib VFPLib)

add_test(NAME parallel-shock COMMAND
  mpirun -np ${MPI_NPROC} ./test-parallel-shock parameter.prm 0.33
)
set_property(TEST parallel-shock PROPERTY LABELS
  "parallel-shock"
  "2d" "momentum" "spatial_advection" "source" "scattering"
  "shock" "CN" "mpi"
  "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
