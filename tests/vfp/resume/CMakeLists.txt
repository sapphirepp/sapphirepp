file(COPY parameter-create-checkpoint.prm
  parameter-restart-new-checkpoint.prm
  parameter-restart-old-checkpoint.prm
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)
file(COPY old-checkpoint DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


add_executable(test-resume-vfp test-resume-vfp.cpp)
sapphirepp_setup_target_vfp(test-resume-vfp
  ${PROJECT_SOURCE_DIR}/examples/vfp/parallel-shock
)
target_include_directories(test-resume-vfp PUBLIC
  ${PROJECT_SOURCE_DIR}/tests/vfp
)

add_test(NAME create-checkpoint-vfp COMMAND
  ./test-resume-vfp parameter-create-checkpoint.prm
)
set_property(TEST create-checkpoint-vfp PROPERTY FIXTURES_SETUP resume-vfp)
set_property(TEST create-checkpoint-vfp PROPERTY LABELS
  "resume" "checkpoint"
  "2d" "momentum" "spatial_advection" "source" "scattering"
  "shock" "CN" "upwind-flux"
  "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)


add_test(NAME restart-new-checkpoint-vfp COMMAND
  mpirun -np ${MPI_NPROC} ./test-resume-vfp parameter-restart-new-checkpoint.prm --resume
)
set_property(TEST restart-new-checkpoint-vfp PROPERTY FIXTURES_REQUIRED resume-vfp)
set_property(TEST restart-new-checkpoint-vfp PROPERTY LABELS
  "resume" "restart"
  "2d" "momentum" "spatial_advection" "source" "scattering"
  "shock" "CN" "upwind-flux" "mpi"
  "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)


add_test(NAME restart-old-checkpoint-vfp COMMAND
  ./test-resume-vfp parameter-restart-old-checkpoint.prm --resume
)
set_property(TEST restart-old-checkpoint-vfp PROPERTY LABELS
  "resume" "restart"
  "2d" "momentum" "spatial_advection" "source" "scattering"
  "shock" "CN" "upwind-flux"
  "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
