file(COPY parameter.prm
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)


add_executable(test-steady-state-parallel-shock
  test-steady-state-parallel-shock.cpp)
sapphirepp_setup_target_vfp(test-steady-state-parallel-shock
  ${PROJECT_SOURCE_DIR}/examples/vfp/steady-state-parallel-shock
)
target_include_directories(test-steady-state-parallel-shock PUBLIC
  ${PROJECT_SOURCE_DIR}/tests/vfp
)

add_test(NAME steady-state-parallel-shock COMMAND
  ./test-steady-state-parallel-shock parameter.prm 0.33
)
set_property(TEST steady-state-parallel-shock PROPERTY LABELS
  "steady-state-parallel-shock"
  "2d" "momentum" "spatial_advection" "source" "scattering"
  "steady_sate" "shock" "upwind-flux"
  "commit-Release" "commit-Debug"
  "full-suite-Release" "full-suite-Debug"
)


add_executable(test-steady-state-parallel-shock-scaled
  test-steady-state-parallel-shock-scaled.cpp)
sapphirepp_setup_target_vfp(test-steady-state-parallel-shock-scaled
  ${PROJECT_SOURCE_DIR}/examples/vfp/steady-state-parallel-shock/scaled
)
target_include_directories(test-steady-state-parallel-shock-scaled PUBLIC
  ${PROJECT_SOURCE_DIR}/tests/vfp
)

add_test(NAME steady-state-parallel-shock-scaled COMMAND
  ./test-steady-state-parallel-shock-scaled parameter.prm 19
)
set_property(TEST steady-state-parallel-shock-scaled PROPERTY LABELS
  "steady-state-parallel-shock-scaled"
  "2d" "momentum" "spatial_advection" "source" "scattering"
  "steady_sate" "scaled_distribution_function" "shock" "upwind-flux"
  "commit-Release" "commit-Debug"
  "full-suite-Release" "full-suite-Debug"
)


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
