add_executable(parallel-shock parallel-shock.cpp)
sapphirepp_setup_target_vfp(parallel-shock)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
