add_executable(linear-hydro-waves linear-hydro-waves.cpp ${MHD_SOURCES})
deal_ii_setup_target(linear-hydro-waves)
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(linear-hydro-waves PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(linear-hydro-waves PUBLIC
                           ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd)
target_link_libraries(linear-hydro-waves UtilsLib)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
