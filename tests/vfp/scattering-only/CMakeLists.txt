file(COPY parameter-1d-CN.prm parameter-1d-ERK4.prm parameter-1d-LSERK4.prm
  parameter-1d-FE.prm parameter-1d-BE.prm
  parameter-2d-ERK4.prm
  parameter-3d-ERK4.prm
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)


add_executable(test-scattering-only-1d test-scattering-only-1d.cpp)
sapphirepp_setup_target_vfp(test-scattering-only-1d
  ${PROJECT_SOURCE_DIR}/examples/vfp/scattering-only
)
target_include_directories(test-scattering-only-1d PUBLIC
  ${PROJECT_SOURCE_DIR}/tests/vfp
)

add_test(NAME scattering-only-1d-CN COMMAND
  ./test-scattering-only-1d parameter-1d-CN.prm 1e-2
)
set_property(TEST scattering-only-1d-CN PROPERTY LABELS
  "scattering-only"
  "1d" "scattering"
  "CN" "periodic BC"
  "full-suite-Release" "full-suite-Debug"
)

add_test(NAME scattering-only-1d-ERK4 COMMAND
  ./test-scattering-only-1d parameter-1d-ERK4.prm 1e-4
)
set_property(TEST scattering-only-1d-ERK4 PROPERTY LABELS
  "scattering-only"
  "1d" "scattering"
  "ERK4" "periodic BC"
  "fast" "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)

add_test(NAME scattering-only-1d-LSERK4 COMMAND
  ./test-scattering-only-1d parameter-1d-LSERK4.prm 3e-3
)
set_property(TEST scattering-only-1d-LSERK4 PROPERTY LABELS
  "scattering-only"
  "1d" "scattering"
  "LSERK4" "periodic BC"
  "full-suite-Debug" "full-suite-Release"
)

add_test(NAME scattering-only-1d-FE COMMAND
  ./test-scattering-only-1d parameter-1d-FE.prm 1e-1
)
set_property(TEST scattering-only-1d-FE PROPERTY LABELS
  "scattering-only"
  "1d" "scattering"
  "FE" "periodic BC"
  "full-suite-Debug" "full-suite-Release"
)

add_test(NAME scattering-only-1d-BE COMMAND
  ./test-scattering-only-1d parameter-1d-BE.prm 1e-1
)
set_property(TEST scattering-only-1d-BE PROPERTY LABELS
  "scattering-only"
  "1d" "scattering"
  "BE" "periodic BC"
  "full-suite-Debug" "full-suite-Release"
)


add_executable(test-scattering-only-2d test-scattering-only-2d.cpp)
sapphirepp_setup_target_vfp(test-scattering-only-2d
  ${PROJECT_SOURCE_DIR}/examples/vfp/scattering-only
)
target_include_directories(test-scattering-only-2d PUBLIC
  ${PROJECT_SOURCE_DIR}/tests/vfp
)

add_test(NAME scattering-only-2d COMMAND
  ./test-scattering-only-2d parameter-2d-ERK4.prm 1e-4
)
set_property(TEST scattering-only-2d PROPERTY LABELS
  "scattering-only"
  "2d" "scattering"
  "ERK4" "periodic BC"
  "full-suite-Release"
)


add_executable(test-scattering-only-3d test-scattering-only-3d.cpp)
sapphirepp_setup_target_vfp(test-scattering-only-3d
  ${PROJECT_SOURCE_DIR}/examples/vfp/scattering-only
)
target_include_directories(test-scattering-only-3d PUBLIC
  ${PROJECT_SOURCE_DIR}/tests/vfp
)

add_test(NAME scattering-only-3d COMMAND
  ./test-scattering-only-3d parameter-3d-ERK4.prm 1e-4
)
set_property(TEST scattering-only-3d PROPERTY LABELS
  "scattering-only"
  "3d" "scattering"
  "ERK4" "periodic BC"
  "commit-Release"
  "full-suite-Release"
)


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
