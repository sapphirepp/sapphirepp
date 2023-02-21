add_library(VFPSolverLib reference-values.cpp physical-setup.cpp particle-functions.cpp pde-system.cpp vfp-solver-control.cpp upwind-flux.cpp)
deal_ii_setup_target(VFPSolverLib)
target_include_directories(VFPSolverLib PUBLIC ${PROJECT_SOURCE_DIR}/include/)
