add_library(MHDLib
  mhd-equations.cpp
  mhd-parameters.cpp
  numerical-flux.cpp
  slope-limiter.cpp
)
deal_ii_setup_target(MHDLib)
target_include_directories(MHDLib PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/utils
)
target_include_directories(MHDLib PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd
)
target_link_libraries(MHDLib UtilsLib)
