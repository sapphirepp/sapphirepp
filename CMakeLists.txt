cmake_minimum_required(VERSION 3.5.0)
# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)

set(DEAL_II_MIN_VERSION "9.4.0") #TODO: check if this is the correct version
set(DOXYGEN_MIN_VERSION "1.8.0") #TODO: check if this is the correct version

file(READ ${CMAKE_SOURCE_DIR}/VERSION PROJECT_VERSION)
message("Sapphire++ version ${PROJECT_VERSION}")

option(SAPPHIRE "Compile Sapphire++" ON )
option(VFP "Build with VFP module" ON)
option(EXAMPLES "Build with examples" OFF)
option(TESTS "Build with tests" OFF)
option(DOC "Build with documentation" OFF)

project(Sapphire VERSION ${PROJECT_VERSION})

if(SAPPHIRE)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # necessary for the language server

  set(DEAL_II_DIR "~/.local/lib/dealii")


  find_package(deal.II ${DEAL_II_MIN_VERSION} REQUIRED
               HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

  deal_ii_initialize_cached_variables()
  file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR})

  configure_file(${PROJECT_SOURCE_DIR}/include/sapphirepp/version.h.in
                 ${PROJECT_BINARY_DIR}/include/sapphirepp/version.h @ONLY)

  add_subdirectory(src/utils)

  if(VFP)
    # VFP source files have to be compiled together with the main/config file
    file(GLOB VFP_SOURCES src/vfp/*.cpp)

    file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/include)
    add_executable(sapphire sapphire.cpp ${VFP_SOURCES})
    deal_ii_setup_target(sapphire)
    target_include_directories(sapphire PUBLIC
                               ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
    target_include_directories(sapphire PUBLIC ${CONFIG_HEADER_DIR})
    target_link_libraries(sapphire UtilsLib)
  endif(VFP)

  if(EXAMPLES)
    add_subdirectory(examples)
  endif(EXAMPLES)

  if(TESTS)
    include(CTest)
    add_subdirectory(tests)
  endif(TESTS)
endif(SAPPHIRE)

  
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
  endif(DOXYGEN_FOUND)
endif(DOC)
