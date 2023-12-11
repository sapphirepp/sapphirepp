if(NOT DEFINED ENV{MPI_NPROC})
  set(MPI_NPROC 8)
else()
  set(MPI_NPROC $ENV{MPI_NPROC})
endif()

add_subdirectory(parallel-shock)
add_subdirectory(gyro-motion-f0)
add_subdirectory(scattering-only)
