add_library(VFPSolverLib reference-values.cpp)
deal_ii_setup_target(VFPSolverLib)
target_include_directories(VFPSolverLib PUBLIC ${PROJECT_SOURCE_DIR}/include/)
