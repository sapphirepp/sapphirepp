# Define the executable
add_executable(scattering-only scattering-only.cpp)

# Automatically link to Sapphire++ and deal.II. 
# The second argument specifies the folder with the config.h file
# in order to include it in the library.
# (Can be omitted if it is equal to ${CMAKE_CURRENT_SOURCE_DIR})
sapphirepp_setup_target_vfp(scattering-only ${CMAKE_CURRENT_SOURCE_DIR})

# Copy the parameter file to the build directory
file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Clean up the results folder when `make clean` is called
set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
