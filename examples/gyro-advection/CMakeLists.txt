add_executable(gyro-advection gyro-advection.cpp ${VFP_SOURCES})
deal_ii_setup_target(gyro-advection)
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(gyro-advection PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(gyro-advection PUBLIC
                           ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
target_link_libraries(gyro-advection UtilsLib)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
