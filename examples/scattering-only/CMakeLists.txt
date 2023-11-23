# Define the executable. Notice that the ${VFP_SOURCES} need to be included.
add_executable(scattering-only scattering-only.cpp ${VFP_SOURCES})

# Link to deal.II
deal_ii_setup_target(scattering-only)

# We have to include the current folder to include the config.h file
file(GLOB CONFIG_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(scattering-only PUBLIC ${CONFIG_HEADER_DIR})

# Link the rest of Sapphire++
target_include_directories(scattering-only PUBLIC 
                           ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
target_link_libraries(scattering-only UtilsLib)

# Copy the parameter file to the build directory
file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Clean up the results folder when `make clean` is called
set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
