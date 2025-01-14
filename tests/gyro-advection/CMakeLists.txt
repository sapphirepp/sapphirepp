file(COPY parameter-advection.prm parameter-static.prm
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)


add_executable(test-gyro-advection test-gyro-advection.cpp ${VFP_SOURCES})
deal_ii_setup_target(test-gyro-advection)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/gyro-advection)
target_include_directories(test-gyro-advection PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp
)
target_include_directories(test-gyro-advection PUBLIC ${CONFIG_HEADER_DIR})
target_link_libraries(test-gyro-advection UtilsLib)

add_test(NAME gyro-advection COMMAND
  ./test-gyro-advection parameter-advection.prm 0.15
)
set_property(TEST gyro-advection PROPERTY LABELS
  "gyro-advection"
  "2d" "spatial_advection" "magnetic"
  "CN" "periodic BC"
  "commit-Release" "commit-Debug"
  "full-suite-Release" "full-suite-Debug"
)

add_test(NAME gyro-static COMMAND
  ./test-gyro-advection parameter-static.prm 5e-3
)
set_property(TEST gyro-static PROPERTY LABELS
  "gyro-advection"
  "2d" "spatial_advection" "magnetic"
  "CN" "continuos gradients"
  "commit-Release" "commit-Debug"
  "full-suite-Release" "full-suite-Debug"
)


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
