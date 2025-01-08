add_library(UtilsLib
  output-parameters.cpp
  sapphirepp-logstream.cpp
  tools.cpp)
deal_ii_setup_target(UtilsLib)
target_include_directories(UtilsLib PUBLIC
                           ${PROJECT_BINARY_DIR}/include/sapphirepp)
target_include_directories(UtilsLib PUBLIC
                           ${PROJECT_SOURCE_DIR}/include/sapphirepp/utils)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost REQUIRED COMPONENTS program_options)
if(Boost_FOUND)
  target_include_directories(UtilsLib PUBLIC ${Boost_INCLUDE_DIRS})
  target_link_libraries(UtilsLib ${Boost_LIBRARIES})
endif()
