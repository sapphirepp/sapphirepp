add_executable(closure closure.cpp)
sapphirepp_setup_target_vfp(closure)

file(COPY parameter.prm heatmap.gp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
