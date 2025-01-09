add_executable(convergence-study convergence-study.cpp ${VFP_SOURCES})
deal_ii_setup_target(convergence-study)
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(convergence-study PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(convergence-study PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp
)
target_link_libraries(convergence-study UtilsLib VFPLib)

file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
