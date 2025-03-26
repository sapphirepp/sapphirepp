add_library(VFPLib
  particle-functions.cpp
  pde-system.cpp
  probe-location.cpp
  vfp-parameters.cpp
)
deal_ii_setup_target(VFPLib)
target_include_directories(VFPLib PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/utils
)
target_include_directories(VFPLib PUBLIC
  ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp
)
target_link_libraries(VFPLib UtilsLib)
