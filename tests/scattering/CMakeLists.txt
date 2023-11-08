file(COPY parameter-1d-CN.prm parameter-1d-ERK4.prm parameter-1d-LSERK.prm
          parameter-1d-FE.prm parameter-1d-BE.prm
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(scattering_test scattering_test.cpp ${VFP_SOURCES})
deal_ii_setup_target(scattering_test)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/scattering)
target_include_directories(scattering_test PUBLIC ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
target_include_directories(scattering_test PUBLIC ${CONFIG_HEADER_DIR})
target_link_libraries(scattering_test UtilsLib)

add_test(NAME scattering-1d-CN
         COMMAND ./scattering_test parameter-1d-CN.prm test_run 2e-2)
set_property(TEST scattering-1d-CN PROPERTY
             LABELS "scattering"
                    "Crank-Nicholson" "1d" "periodic BC")

add_test(NAME scattering-1d-ERK4
         COMMAND ./scattering_test parameter-1d-ERK4.prm test_run 2e-4)
set_property(TEST scattering-1d-ERK4 PROPERTY
             LABELS "scattering"
                    "ERK4" "1d" "periodic BC"
                    "fast")

# TODO: LSERK error larger than expected
add_test(NAME scattering-1d-LSERK
         COMMAND ./scattering_test parameter-1d-LSERK.prm test_run 2e-2)
set_property(TEST scattering-1d-LSERK PROPERTY
             LABELS "scattering"
                    "LSERK" "1d"  "periodic BC")

add_test(NAME scattering-1d-FE
         COMMAND ./scattering_test parameter-1d-FE.prm test_run 3e-1)
set_property(TEST scattering-1d-FE PROPERTY
             LABELS "scattering"
                    "Forward Euler" "1d" "periodic BC")

add_test(NAME scattering-1d-BE
         COMMAND ./scattering_test parameter-1d-BE.prm test_run 3e-1)
set_property(TEST scattering-1d-BE PROPERTY
             LABELS "scattering"
                    "Backward Euler" "1d" "periodic BC")


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
