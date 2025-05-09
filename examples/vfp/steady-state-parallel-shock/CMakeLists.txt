add_executable(steady-state-parallel-shock steady-state-parallel-shock.cpp)
sapphirepp_setup_target_vfp(steady-state-parallel-shock)

add_subdirectory(scaled)
