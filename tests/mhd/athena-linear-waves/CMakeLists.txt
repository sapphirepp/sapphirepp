file(COPY parameter-1d-sound.prm parameter-1d-entropy.prm
  parameter-1d-fast.prm parameter-1d-alfven.prm parameter-1d-slow.prm
  parameter-2d-sound.prm parameter-2d-alfven.prm
  parameter-3d-sound.prm parameter-3d-alfven.prm
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)


add_executable(test-athena-linear-waves-1d test-athena-linear-waves-1d.cpp ${MHD_SOURCES})
deal_ii_setup_target(test-athena-linear-waves-1d)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/mhd/athena-linear-waves)
target_include_directories(test-athena-linear-waves-1d PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(test-athena-linear-waves-1d PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd
)
target_link_libraries(test-athena-linear-waves-1d UtilsLib MHDLib)

add_test(NAME athena-linear-waves-1d-sound COMMAND
  ./test-athena-linear-waves-1d parameter-1d-sound.prm 4e-7
)
set_property(TEST athena-linear-waves-1d-sound PROPERTY LABELS
  "athena-linear-waves"
  "1d" "mhd" "sound-wave"
  "ERK4" "periodic BC"
  "fast" "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)

add_test(NAME athena-linear-waves-1d-entropy COMMAND
  ./test-athena-linear-waves-1d parameter-1d-entropy.prm 2e-7
)
set_property(TEST athena-linear-waves-1d-entropy PROPERTY LABELS
  "athena-linear-waves"
  "1d" "mhd" "entropy-wave"
  "ERK4" "periodic BC"
  "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)

add_test(NAME athena-linear-waves-1d-fast COMMAND
  ./test-athena-linear-waves-1d parameter-1d-fast.prm 2e-7
)
set_property(TEST athena-linear-waves-1d-fast PROPERTY LABELS
  "athena-linear-waves"
  "1d" "mhd" "fast-wave"
  "ERK4" "periodic BC"
  "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)

add_test(NAME athena-linear-waves-1d-alfven COMMAND
  ./test-athena-linear-waves-1d parameter-1d-alfven.prm 8e-8
)
set_property(TEST athena-linear-waves-1d-alfven PROPERTY LABELS
  "athena-linear-waves"
  "1d" "mhd" "alfven-wave"
  "ERK4" "periodic BC"
  "fast" "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)

add_test(NAME athena-linear-waves-1d-slow COMMAND
  ./test-athena-linear-waves-1d parameter-1d-slow.prm 8e-8
)
set_property(TEST athena-linear-waves-1d-slow PROPERTY LABELS
  "athena-linear-waves"
  "1d" "mhd" "slow-wave"
  "ERK4" "periodic BC"
  "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)


add_executable(test-athena-linear-waves-2d test-athena-linear-waves-2d.cpp ${MHD_SOURCES})
deal_ii_setup_target(test-athena-linear-waves-2d)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/mhd/athena-linear-waves)
target_include_directories(test-athena-linear-waves-2d PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(test-athena-linear-waves-2d PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd
)
target_link_libraries(test-athena-linear-waves-2d UtilsLib MHDLib)

add_test(NAME athena-linear-waves-2d-sound COMMAND
  ./test-athena-linear-waves-2d parameter-2d-sound.prm 2e-6
)
set_property(TEST athena-linear-waves-2d-sound PROPERTY LABELS
  "athena-linear-waves"
  "2d" "mhd" "sound-wave"
  "ERK4" "periodic BC"
  "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)

add_test(NAME athena-linear-waves-2d-alfven COMMAND
  ./test-athena-linear-waves-2d parameter-2d-alfven.prm 5e-7
)
set_property(TEST athena-linear-waves-2d-alfven PROPERTY LABELS
  "athena-linear-waves"
  "2d" "mhd" "alfven-wave"
  "ERK4" "periodic BC"
  "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)


add_executable(test-athena-linear-waves-3d test-athena-linear-waves-3d.cpp ${MHD_SOURCES})
deal_ii_setup_target(test-athena-linear-waves-3d)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/mhd/athena-linear-waves)
target_include_directories(test-athena-linear-waves-3d PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(test-athena-linear-waves-3d PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd
)
target_link_libraries(test-athena-linear-waves-3d UtilsLib MHDLib)

add_test(NAME athena-linear-waves-3d-sound COMMAND
  ./test-athena-linear-waves-3d parameter-3d-sound.prm 2e-5
)
set_property(TEST athena-linear-waves-3d-sound PROPERTY LABELS
  "athena-linear-waves"
  "3d" "mhd" "sound-wave"
  "ERK4" "periodic BC"
  "full-suite-Release"
)

add_test(NAME athena-linear-waves-3d-alfven COMMAND
  ./test-athena-linear-waves-3d parameter-3d-alfven.prm 6e-6
)
set_property(TEST athena-linear-waves-3d-alfven PROPERTY LABELS
  "athena-linear-waves"
  "3d" "mhd" "alfven-wave"
  "ERK4" "periodic BC"
  "commit-Release"
  "full-suite-Release"
)


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
