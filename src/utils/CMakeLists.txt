add_library(SapphireUtilsLib
  parameter_parser.cpp)
deal_ii_setup_target(SapphireUtilsLib)
target_include_directories(SapphireUtilsLib PUBLIC ${PROJECT_SOURCE_DIR}/include/utils)
