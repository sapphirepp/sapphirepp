add_executable(parallel-shock parallel-shock.cpp ${VFP_SOURCES})
deal_ii_setup_target(parallel-shock)
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(parallel-shock PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(parallel-shock PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp
)
target_link_libraries(parallel-shock UtilsLib VFPLib)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
