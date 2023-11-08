add_executable(scattering-only scattering-only.cpp ${VFP_SOURCES})
deal_ii_setup_target(scattering-only)
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(scattering-only PUBLIC ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
target_include_directories(scattering-only PUBLIC ${CONFIG_HEADER_DIR})
target_link_libraries(scattering-only UtilsLib)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
