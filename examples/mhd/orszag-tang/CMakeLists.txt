add_executable(orszag-tang orszag-tang.cpp ${MHD_SOURCES})
deal_ii_setup_target(orszag-tang)
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(orszag-tang PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(orszag-tang PUBLIC
                           ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd)
target_link_libraries(orszag-tang UtilsLib)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
