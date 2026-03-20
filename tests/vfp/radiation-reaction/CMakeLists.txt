file(COPY parameter.prm parameter-massratio.prm
  parameter-with-rotation.prm parameter-with-source.prm
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)


add_executable(test-radiation-reaction test-radiation-reaction.cpp)
sapphirepp_setup_target_vfp(test-radiation-reaction
  ${PROJECT_SOURCE_DIR}/examples/vfp/radiation-reaction/radiation-reaction
)
target_include_directories(test-radiation-reaction PUBLIC
  ${PROJECT_SOURCE_DIR}/tests/vfp
)

add_test(NAME radiation-reaction COMMAND
  ./test-radiation-reaction parameter.prm 5e-6
)
set_property(TEST radiation-reaction PROPERTY LABELS
  "radiation-reaction"
  "1d" "momentum" "spatial_advection" "radiation_reaction"
  "scaled_distribution_function"
  "ERK4" "zero-inflow" "upwind-flux"
  "fast" "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)

add_test(NAME radiation-reaction-massratio COMMAND
  ./test-radiation-reaction parameter-massratio.prm 5e-6
)
set_property(TEST radiation-reaction-massratio PROPERTY LABELS
  "radiation-reaction-massratio"
  "1d" "momentum" "spatial_advection" "radiation_reaction"
  "scaled_distribution_function"
  "ERK4" "zero-inflow" "upwind-flux"
  "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)


add_executable(test-radiation-reaction-with-rotation test-radiation-reaction-with-rotation.cpp)
sapphirepp_setup_target_vfp(test-radiation-reaction-with-rotation
  ${PROJECT_SOURCE_DIR}/examples/vfp/radiation-reaction/radiation-reaction-with-rotation
)
target_include_directories(test-radiation-reaction-with-rotation PUBLIC
  ${PROJECT_SOURCE_DIR}/tests/vfp
)

add_test(NAME radiation-reaction-with-rotation COMMAND
  ./test-radiation-reaction-with-rotation parameter-with-rotation.prm 0.01
)
set_property(TEST radiation-reaction-with-rotation PROPERTY LABELS
  "radiation-reaction-with-rotation"
  "1d" "momentum" "spatial_advection" "radiation_reaction" "rotation"
  "scaled_distribution_function"
  "CN" "zero-inflow" "upwind-flux"
  "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)

add_executable(test-radiation-reaction-with-source test-radiation-reaction-with-source.cpp)
sapphirepp_setup_target_vfp(test-radiation-reaction-with-source
  ${PROJECT_SOURCE_DIR}/examples/vfp/radiation-reaction/radiation-reaction-with-source
)
target_include_directories(test-radiation-reaction-with-source PUBLIC
  ${PROJECT_SOURCE_DIR}/tests/vfp
)

add_test(NAME radiation-reaction-with-source COMMAND
  ./test-radiation-reaction-with-source parameter-with-source.prm 0.01
)
set_property(TEST radiation-reaction-with-source PROPERTY LABELS
  "radiation-reaction-with-source"
  "1d" "momentum" "spatial_advection" "radiation_reaction" "source"
  "scaled_distribution_function"
  "LSERK4" "zero-inflow" "local_lax_friedrichs_flux"
  "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
