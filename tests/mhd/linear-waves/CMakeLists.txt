file(COPY parameter-1d-alfven.prm 
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test-linear-waves-1d test-linear-waves-1d.cpp ${MHD_SOURCES})
deal_ii_setup_target(test-linear-waves-1d)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/mhd/linear-waves)
target_include_directories(test-linear-waves-1d PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(test-linear-waves-1d PUBLIC
                           ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd)
target_link_libraries(test-linear-waves-1d UtilsLib)

add_test(NAME linear-waves-1d-alfven
         COMMAND ./test-linear-waves-1d parameter-1d-alfven.prm 8e-8)
set_property(TEST linear-waves-1d-alfven PROPERTY
             LABELS "linear-waves"
                    "1d" "mhd" "alfven-wave"
                    "ERK4" "periodic BC"
                    "fast" "commit-Debug" "commit-Release"
                    "full-suite-Release" "full-suite-Debug")

                    
set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
