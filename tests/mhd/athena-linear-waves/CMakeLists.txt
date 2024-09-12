file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test-athena-linear-waves test-athena-linear-waves.cpp ${MHD_SOURCES})
deal_ii_setup_target(test-athena-linear-waves)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/mhd/athena-linear-waves)
target_include_directories(test-athena-linear-waves PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(test-athena-linear-waves PUBLIC
                           ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd)
target_link_libraries(test-athena-linear-waves UtilsLib)


add_test(NAME athena-linear-waves
         COMMAND ./test-athena-linear-waves parameter.prm 1e-5)
set_property(TEST athena-linear-waves PROPERTY
             LABELS "athena-linear-waves"
                    "1d" "mhd"
                    "ERK4" "periodic BC"
                    "full-suite-Release" "full-suite-Debug")


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
