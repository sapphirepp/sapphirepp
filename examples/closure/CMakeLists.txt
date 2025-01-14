add_executable(closure closure.cpp ${VFP_SOURCES})
deal_ii_setup_target(closure)
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(closure PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(closure PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp
)
target_link_libraries(closure UtilsLib)

file(COPY parameter.prm heatmap.gp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
