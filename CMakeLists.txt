cmake_minimum_required(VERSION 3.1.0)
# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)

# if( NOT CMAKE_BUILD_TYPE )       # be explicit, if you want a debug build
#   set( CMAKE_BUILD_TYPE Debug) # CMAKE_BUILD_TYPE information is also
# endif()   
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # necessary for the language server (e.g. to recognise the include libraries)

set(DEAL_II_DIR "~/.local/lib/dealii")

find_package(deal.II 9.4.0 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
deal_ii_initialize_cached_variables() # CMAKE_BUILD_TYPE is set to the whatever dealII was built with

project(Sapphire VERSION 1.0.2)
add_subdirectory(src/utils)
option(Sapphire "Compile Sapphire++" ON )
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR})
if(Sapphire)
  file(GLOB VFP_SOURCES src/vfp/*.cpp)
  add_executable(sapphire sapphire.cpp ${VFP_SOURCES})
  deal_ii_setup_target(sapphire)
  target_include_directories(sapphire PUBLIC ${PROJECT_SOURCE_DIR}/include/vfp)
  target_include_directories(sapphire PUBLIC ${CONFIG_HEADER_DIR})
  target_link_libraries(sapphire UtilsLib)
endif(Sapphire)

find_package(Doxygen)
if(DOXYGEN_FOUND)
    # Configuration file for Doxygen
    set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/doc/Doxyfile.in)
    set(DOXYGEN_OUT ${PROJECT_BINARY_DIR}/doc/Doxyfile)

    # Generate the Doxyfile from the template
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    # Add custom target for Doxygen documentation
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Generating Doxygen documentation"
        VERBATIM
    )
endif()


