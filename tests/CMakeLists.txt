if(NOT DEFINED ENV{MPI_NPROC})
  set(MPI_NPROC 8)
else()
  set(MPI_NPROC $ENV{MPI_NPROC})
endif()

add_test(NAME quick-start
         COMMAND ../sapphirepp)
set_property(TEST quick-start PROPERTY
             LABELS "quick-start"
                    "2d" "momentum" "spatial_advection" "source" "scattering"
                    "shock" "Crank-Nicholson"
                    "commit-Release" "commit-Debug"
                    "full-suite-Release" "full-suite-Debug")

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
"${CMAKE_CURRENT_BINARY_DIR}/results"
)

add_subdirectory(parallel-shock)
add_subdirectory(gyro-advection)
add_subdirectory(scattering-only)
