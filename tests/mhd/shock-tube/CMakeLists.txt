file(COPY parameter-sod.prm
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)


add_executable(test-shock-tube test-shock-tube.cpp)
sapphirepp_setup_target_mhd(test-shock-tube
  ${PROJECT_SOURCE_DIR}/examples/mhd/shock-tube
)
target_include_directories(test-shock-tube PUBLIC
  ${PROJECT_SOURCE_DIR}/tests/mhd
)

add_test(NAME shock-tube-sod COMMAND
  ./test-shock-tube parameter-sod.prm
  # "athena-results"
  "/home/schulze/Documents/PhD/Code/athena-results" # TODO: copy data over
  "Sod"
  6e-3
)
set_property(TEST shock-tube-sod PROPERTY LABELS
  "shock-tube"
  "1d" "sod"
  "ERK2" "zero-inflow"
  "fast" "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
