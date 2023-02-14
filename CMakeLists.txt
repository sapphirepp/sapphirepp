cmake_minimum_required(VERSION 3.1.0)
# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)

# if( NOT CMAKE_BUILD_TYPE )       # be explicit, if you want a debug build
#   set( CMAKE_BUILD_TYPE Debug) # CMAKE_BUILD_TYPE information is also
# endif()   
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # necessary for the language server (e.g. to recognise the include libraries)

set(DEAL_II_DIR "~/.local/lib/dealii")

find_package(deal.II 9.3.0 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
deal_ii_initialize_cached_variables() # CMAKE_BUILD_TYPE is set to the whatever dealII was built with

project(VFP-Equation)
option(vfp-equation-solver "Compile the VFP Equation Solver" ON )
if(vfp-equation-solver)
  add_subdirectory(src)
  add_executable(vfp-equation-solver vfp-equation-solver.cpp)
  deal_ii_setup_target(vfp-equation-solver)
  target_link_libraries(vfp-equation-solver VFPSolverLib)
endif(vfp-equation-solver)



