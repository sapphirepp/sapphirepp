add_library(SapphireUtilsLib
  output_module.cpp)
deal_ii_setup_target(SapphireUtilsLib)
target_include_directories(SapphireUtilsLib PUBLIC ${PROJECT_SOURCE_DIR}/include/utils)
