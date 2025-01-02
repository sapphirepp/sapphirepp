add_executable(steady-state-parallel-shock steady-state-parallel-shock.cpp ${VFP_SOURCES})
deal_ii_setup_target(steady-state-parallel-shock)
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(steady-state-parallel-shock PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(steady-state-parallel-shock PUBLIC
                           ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
target_link_libraries(steady-state-parallel-shock UtilsLib)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
