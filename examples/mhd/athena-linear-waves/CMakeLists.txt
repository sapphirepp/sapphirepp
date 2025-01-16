add_executable(athena-linear-waves athena-linear-waves.cpp)
sapphirepp_setup_target_mhd(athena-linear-waves)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
