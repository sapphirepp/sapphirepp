add_executable(shock-tube-2d shock-tube-2d.cpp)
sapphirepp_setup_target_mhd(shock-tube-2d)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
