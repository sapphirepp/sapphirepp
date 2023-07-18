add_library(SapphireUtilsLib
  output-module.cpp
  parameter-parser.cpp)
deal_ii_setup_target(SapphireUtilsLib)
target_include_directories(SapphireUtilsLib PUBLIC ${PROJECT_SOURCE_DIR}/include/utils)
