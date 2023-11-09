add_executable(gyro-motion-f0 gyro-motion-f0.cpp ${VFP_SOURCES})
deal_ii_setup_target(gyro-motion-f0)
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(gyro-motion-f0 PUBLIC ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
target_include_directories(gyro-motion-f0 PUBLIC ${CONFIG_HEADER_DIR})
target_link_libraries(gyro-motion-f0 UtilsLib)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
