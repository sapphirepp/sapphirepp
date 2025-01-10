add_executable(linear-waves linear-waves.cpp ${MHD_SOURCES})
deal_ii_setup_target(linear-waves)
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(linear-waves PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(linear-waves PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd
)
target_link_libraries(linear-waves UtilsLib MHDLib)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
