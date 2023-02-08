CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)
# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# if( NOT CMAKE_BUILD_TYPE )       # be explicit, if you want a debug build
#   set( CMAKE_BUILD_TYPE Debug) # CMAKE_BUILD_TYPE information is also
# endif()   
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # necessary for the language server (e.g. to recognise the include libraries)

SET(DEAL_II_DIR "~/.local/lib/dealii")

FIND_PACKAGE(deal.II 9.3.0 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
DEAL_II_INITIALIZE_CACHED_VARIABLES() # CMAKE_BUILD_TYPE is set to the whatever dealII was built with

PROJECT(VFP-Equation)
option(vfp-equation-solver "Compile the pure advection system" ON )
if(vfp-equation-solver)
  ADD_EXECUTABLE(vfp-equation-solver vfp-equation-solver.cpp)
  DEAL_II_SETUP_TARGET(vfp-equation-solver)
endif(vfp-equation-solver)



