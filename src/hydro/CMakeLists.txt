add_library(HydroSolverLib
  burgers-eq.cpp
  conservation-eq.cpp
  numerics.cpp)
deal_ii_setup_target(HydroSolverLib)
target_include_directories(HydroSolverLib PUBLIC ${PROJECT_SOURCE_DIR}/include/hydro)
target_link_libraries(HydroSolverLib SapphireUtilsLib)
