# Define the executable
add_executable(non-homogeneous-bc non-homogeneous-bc.cpp)

# Automatically link to Sapphire++ and deal.II. 
# The second argument specifies the folder with the config.h file
# in order to include it in the library.
# (Can be omitted if it is equal to ${CMAKE_CURRENT_SOURCE_DIR})
sapphirepp_setup_target_vfp(non-homogeneous-bc ${CMAKE_CURRENT_SOURCE_DIR})
