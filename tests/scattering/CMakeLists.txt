add_executable(scattering_test scattering_test.cpp ${VFP_SOURCES})
deal_ii_setup_target(scattering_test)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/scattering)
target_include_directories(scattering_test PUBLIC ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
target_include_directories(scattering_test PUBLIC ${CONFIG_HEADER_DIR})
target_link_libraries(scattering_test UtilsLib)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME scattering COMMAND ./scattering_test parameter.prm 1e-3)
set_property(TEST scattering PROPERTY LABELS "scattering" "fast" "1d")


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
