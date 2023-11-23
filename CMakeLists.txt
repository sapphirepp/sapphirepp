cmake_minimum_required(VERSION 3.13.0)
# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)

set(DEAL_II_MIN_VERSION "9.5.0")
set(DOXYGEN_MIN_VERSION "1.9.0")

# The version is redefined on the webserver to git describe --tags, where as the
# version file is used for release.
if(NOT DEFINED SAPPHIREPP_VERSION)
  file(READ ${CMAKE_SOURCE_DIR}/VERSION SAPPHIREPP_VERSION)
endif()
string(REGEX REPLACE "([0-9]+\\.[0-9]+\\.[0-9]+)(-.*)?" "\\1"
       PROJECT_VERSION ${SAPPHIREPP_VERSION})

message("Sapphire++ version ${SAPPHIREPP_VERSION}")

option(SAPPHIREPP "Compile Sapphire++" ON )
option(VFP "Build with VFP module" ON)
option(EXAMPLES "Build with examples" OFF)
option(TESTS "Build with tests" OFF)
option(DOC "Build with documentation" OFF)

project(Sapphire++ VERSION ${PROJECT_VERSION})

if(SAPPHIREPP)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # necessary for the language server

  set(DEAL_II_DIR "~/.local/lib/dealii")
  find_package(deal.II ${DEAL_II_MIN_VERSION} REQUIRED
               HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})
  IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake"
    " or set an environment variable \"DEAL_II_DIR\".")
  ENDIF()

  IF(NOT DEAL_II_WITH_PETSC OR NOT DEAL_II_WITH_P4EST)
  MESSAGE(FATAL_ERROR
    "*** deal.II dependencies not met ***\n"
    "Sapphire++ requires deal.II to be configured with the following options:
     DEAL_II_WITH_PETSC = ON
     DEAL_II_WITH_P4EST = ON\n"
    "However, the detected deal.II configuration found at ${DEAL_II_PATH} is:
     DEAL_II_WITH_PETSC = ${DEAL_II_WITH_PETSC}
     DEAL_II_WITH_P4EST = ${DEAL_II_WITH_P4EST}")
  ENDIF()

  deal_ii_initialize_cached_variables()
  file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR})

  configure_file(${PROJECT_SOURCE_DIR}/include/sapphirepp/version.h.in
                 ${PROJECT_BINARY_DIR}/include/sapphirepp/version.h @ONLY)

  add_subdirectory(src/utils)

  if(VFP)
    # VFP source files have to be compiled together with the main/config file
    file(GLOB VFP_SOURCES src/vfp/*.cpp)

    file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/include)
    add_executable(sapphirepp sapphirepp.cpp ${VFP_SOURCES})
    deal_ii_setup_target(sapphirepp)
    target_include_directories(sapphirepp PUBLIC
                               ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
    target_include_directories(sapphirepp PUBLIC ${CONFIG_HEADER_DIR})
    target_link_libraries(sapphirepp UtilsLib)
  endif(VFP)

  if(EXAMPLES)
    add_subdirectory(examples)
  endif(EXAMPLES)

  if(TESTS)
    include(CTest)
    add_subdirectory(tests)
  endif(TESTS)
endif(SAPPHIREPP)

  
if(DOC)
  find_package(Doxygen ${DOXYGEN_MIN_VERSION})
  if(DOXYGEN_FOUND)
      # An empty string implies the usage of a Content Delivery Network (CDN)
      # for Mathjax. Doxygen defaults to https://cdn.jsdelivr.net/npm/mathjax@3
      set(MATHJAX_PATH "" CACHE STRING "MathJax path")

      # Configuration file for Doxygen
      set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/doc/Doxyfile.in)
      set(DOXYGEN_OUT ${PROJECT_BINARY_DIR}/doc/Doxyfile)

      # Generate the Doxyfile from the template
      configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

      # Add custom target for Doxygen documentation
      add_custom_target(doc
          COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
          WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
          COMMENT "Generating Doxygen documentation"
          VERBATIM
      )

    set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
    "${CMAKE_CURRENT_BINARY_DIR}/doc/latest"
    )
  endif(DOXYGEN_FOUND)
endif(DOC)
