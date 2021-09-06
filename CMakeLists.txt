CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)
# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# if( NOT CMAKE_BUILD_TYPE )       # be explicit, if you want a debug build
#   set( CMAKE_BUILD_TYPE Debug) # CMAKE_BUILD_TYPE information is also
# endif()   
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # necessary for the language server (e.g. to recognise the include libraries)

SET(DEAL_II_DIR "/home/schween/.local/lib/dealii")

FIND_PACKAGE(deal.II 9.3.0 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
DEAL_II_INITIALIZE_CACHED_VARIABLES() # CMAKE_BUILD_TYPE is set to the whatever dealII was built with

PROJECT(VFP-Equation)
option(pure_advection_system "Compile the pure advection system" ON )
if(pure_advection_system)
  ADD_EXECUTABLE(pure_advection_system pure_advection_system.cpp)
  DEAL_II_SETUP_TARGET(pure_advection_system)
endif(pure_advection_system)



