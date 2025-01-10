add_executable(shock-tube shock-tube.cpp ${MHD_SOURCES})
deal_ii_setup_target(shock-tube)
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(shock-tube PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(shock-tube PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd
)
target_link_libraries(shock-tube UtilsLib MHDLib)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
