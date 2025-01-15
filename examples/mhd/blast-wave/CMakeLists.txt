add_executable(blast-wave blast-wave.cpp ${MHD_SOURCES})
deal_ii_setup_target(blast-wave)
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(blast-wave PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(blast-wave PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd
)
target_link_libraries(blast-wave UtilsLib)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
