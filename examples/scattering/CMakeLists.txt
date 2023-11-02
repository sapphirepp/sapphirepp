add_executable(scattering scattering.cpp ${VFP_SOURCES})
deal_ii_setup_target(scattering)
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(scattering PUBLIC ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
target_include_directories(scattering PUBLIC ${CONFIG_HEADER_DIR})
target_link_libraries(scattering UtilsLib)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
