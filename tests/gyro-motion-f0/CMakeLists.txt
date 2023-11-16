file(COPY parameter.prm
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


add_executable(test-gyro-motion-f0 test-gyro-motion-f0.cpp ${VFP_SOURCES})
deal_ii_setup_target(test-gyro-motion-f0)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/gyro-motion-f0)
target_include_directories(test-gyro-motion-f0 PUBLIC ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
target_include_directories(test-gyro-motion-f0 PUBLIC ${CONFIG_HEADER_DIR})
target_link_libraries(test-gyro-motion-f0 UtilsLib)

add_test(NAME gyro-motion-f0
         COMMAND ./test-gyro-motion-f0 parameter.prm 3e-2)
set_property(TEST gyro-motion-f0 PROPERTY
             LABELS "gyro-motion-f0"
                    "2d" "spatial_advection" "magnetic"
                    "Crank-Nicholson" "continuos gradients"
                    "commit")


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
