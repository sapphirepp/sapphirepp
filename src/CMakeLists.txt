add_library(VFPSolverLib reference-values.cpp)
target_include_directories(VFPSolverLib PRIVATE ${PROJECT_SOURCE_DIR}/include)
DEAL_II_SETUP_TARGET(VFPSolverLib)
