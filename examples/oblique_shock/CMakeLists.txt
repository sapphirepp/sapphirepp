add_executable(oblique_shock oblique_shock.cpp ${VFP_SOURCES})
deal_ii_setup_target(oblique_shock)
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(oblique_shock PUBLIC ${PROJECT_SOURCE_DIR}/include/vfp)
target_include_directories(oblique_shock PUBLIC ${CONFIG_HEADER_DIR})
target_link_libraries(oblique_shock UtilsLib)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
