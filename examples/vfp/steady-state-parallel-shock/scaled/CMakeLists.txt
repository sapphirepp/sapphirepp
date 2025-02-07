add_executable(steady-state-parallel-shock-scaled steady-state-parallel-shock.cpp)
sapphirepp_setup_target_vfp(steady-state-parallel-shock-scaled)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
