add_library(UtilsLib
  output-parameters.cpp
  sapphirepp-logstream.cpp
  grid-data-function.cpp)
deal_ii_setup_target(UtilsLib)
target_include_directories(UtilsLib PUBLIC
                           ${PROJECT_BINARY_DIR}/include/sapphirepp)
target_include_directories(UtilsLib PUBLIC
                           ${PROJECT_SOURCE_DIR}/include/sapphirepp/utils)
