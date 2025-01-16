add_executable(orszag-tang orszag-tang.cpp)
sapphirepp_setup_target_mhd(orszag-tang)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
