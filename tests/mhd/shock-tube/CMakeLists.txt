file(COPY parameter-sod.prm
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)


add_executable(test-shock-tube test-shock-tube.cpp ${MHD_SOURCES})
deal_ii_setup_target(test-shock-tube)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/mhd/shock-tube)
target_include_directories(test-shock-tube PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(test-shock-tube PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd
  ${PROJECT_SOURCE_DIR}/tests/mhd
)
target_link_libraries(test-shock-tube UtilsLib MHDLib)

add_test(NAME shock-tube-sod COMMAND
  ./test-shock-tube parameter-sod.prm 4e-7
)
set_property(TEST shock-tube-sod PROPERTY LABELS
  "shock-tube"
  "1d" "mhd" "sod"
  "ERK4" "periodic BC"
  "fast" "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
