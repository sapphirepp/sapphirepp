if(NOT DEFINED ENV{MPI_NPROC})
  set(MPI_NPROC 4)
else()
  set(MPI_NPROC $ENV{MPI_NPROC})
endif()

if(SAPPHIREPP_ALL_MODULES)
  add_test(NAME quick-start COMMAND
    ../sapphirepp
  )
  set_property(TEST quick-start PROPERTY LABELS
    "quick-start"
    "vfp"
    "2d" "momentum" "spatial_advection" "source" "scattering"
    "shock" "CN"
    "commit-Release" "commit-Debug"
    "full-suite-Release" "full-suite-Debug"
  )
endif()

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)

if(VFP)
  add_subdirectory(vfp)
endif()

if(MHD)
  add_subdirectory(mhd)
endif()
