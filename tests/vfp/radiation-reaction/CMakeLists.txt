file(COPY parameter-without-rotation.prm
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)


add_executable(test-radiation-reaction-without-rotation test-radiation-reaction-without-rotation.cpp)
sapphirepp_setup_target_vfp(test-radiation-reaction-without-rotation
  ${PROJECT_SOURCE_DIR}/examples/vfp/radiation-reaction
)
target_include_directories(test-radiation-reaction-without-rotation PUBLIC
  ${PROJECT_SOURCE_DIR}/tests/vfp
)

add_test(NAME radiation-reaction-without-rotation COMMAND
  ./test-radiation-reaction-without-rotation parameter-without-rotation.prm 0 # TODO
)
set_property(TEST radiation-reaction-without-rotation PROPERTY LABELS
  "radiation-reaction-without-rotation"
  "1d" "momentum" "spatial_advection" "radiation_reaction"
  "scaled_distribution_function"
  "ERK4" "zero-inflow" "upwind-flux"
  # "fast" "commit-Debug" "commit-Release" # TODO
  "full-suite-Release" "full-suite-Debug"
)


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
