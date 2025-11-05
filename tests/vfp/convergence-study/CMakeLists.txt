file(COPY parameter-CN.prm parameter-ERK4.prm parameter-LSERK4.prm
  parameter-FE.prm parameter-BE.prm
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)


add_executable(test-convergence-study test-convergence-study.cpp)
sapphirepp_setup_target_vfp(test-convergence-study
  ${PROJECT_SOURCE_DIR}/examples/vfp/convergence-study
)
target_include_directories(test-convergence-study PUBLIC
  ${PROJECT_SOURCE_DIR}/tests/vfp
)

add_test(NAME convergence-study-CN COMMAND
  ./test-convergence-study parameter-CN.prm 5e-6
)
set_property(TEST convergence-study-CN PROPERTY LABELS
  "convergence-study"
  "1d" "spatial_advection" "magnetic"
  "CN" "periodic BC" "lax-friedrichs-flux"
  "fast" "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)

add_test(NAME convergence-study-ERK4 COMMAND
  ./test-convergence-study parameter-ERK4.prm 3e-8
)
set_property(TEST convergence-study-ERK4 PROPERTY LABELS
  "convergence-study"
  "1d" "scattering"
  "ERK4" "periodic BC" "lax-friedrichs-flux"
  "fast" "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)

add_test(NAME convergence-study-LSERK4 COMMAND
  ./test-convergence-study parameter-LSERK4.prm 3e-8
)
set_property(TEST convergence-study-LSERK4 PROPERTY LABELS
  "convergence-study"
  "1d" "scattering"
  "LSERK4" "periodic BC" "lax-friedrichs-flux"
  "full-suite-Debug" "full-suite-Release"
)

add_test(NAME convergence-study-FE COMMAND
  ./test-convergence-study parameter-FE.prm 2e-3
)
set_property(TEST convergence-study-FE PROPERTY LABELS
  "convergence-study"
  "1d" "scattering"
  "FE" "periodic BC" "lax-friedrichs-flux"
  "full-suite-Debug" "full-suite-Release"
)

add_test(NAME convergence-study-BE COMMAND
  ./test-convergence-study parameter-BE.prm 6e-4
)
set_property(TEST convergence-study-BE PROPERTY LABELS
  "convergence-study"
  "1d" "scattering"
  "BE" "periodic BC" "lax-friedrichs-flux"
  "full-suite-Debug" "full-suite-Release"
)


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
