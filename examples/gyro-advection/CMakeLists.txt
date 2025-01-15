add_executable(gyro-advection gyro-advection.cpp)
sapphirepp_setup_target_vfp(gyro-advection)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
