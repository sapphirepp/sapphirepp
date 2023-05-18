add_library(HydroSolverLib
  conservation-eq.cpp)
deal_ii_setup_target(HydroSolverLib)
target_include_directories(HydroSolverLib PUBLIC ${PROJECT_SOURCE_DIR}/include/hydro)
