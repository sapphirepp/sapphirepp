file(COPY parameter.prm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test-athena-linear-waves-1d test-athena-linear-waves-1d.cpp ${MHD_SOURCES})
deal_ii_setup_target(test-athena-linear-waves-1d)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/mhd/athena-linear-waves)
target_include_directories(test-athena-linear-waves-1d PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(test-athena-linear-waves-1d PUBLIC
                           ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd)
target_link_libraries(test-athena-linear-waves-1d UtilsLib)

add_test(NAME athena-linear-waves-1d-sound-wave
         COMMAND ./test-athena-linear-waves-1d parameter.prm 1e-5)
set_property(TEST athena-linear-waves-1d-sound-wave PROPERTY
             LABELS "athena-linear-waves"
                    "1d" "mhd" "sound-wave"
                    "ERK4" "periodic BC"
                    "full-suite-Release" "full-suite-Debug")


add_executable(test-athena-linear-waves-2d test-athena-linear-waves-2d.cpp ${MHD_SOURCES})
deal_ii_setup_target(test-athena-linear-waves-2d)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/mhd/athena-linear-waves)
target_include_directories(test-athena-linear-waves-2d PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(test-athena-linear-waves-2d PUBLIC
                           ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd)
target_link_libraries(test-athena-linear-waves-2d UtilsLib)

add_test(NAME athena-linear-waves-2d-sound-wave
         COMMAND ./test-athena-linear-waves-2d parameter.prm 1e-5)
set_property(TEST athena-linear-waves-1d-sound-wave PROPERTY
             LABELS "athena-linear-waves"
                    "2d" "mhd" "sound-wave"
                    "ERK4" "periodic BC"
                    "full-suite-Release" "full-suite-Debug")


add_executable(test-athena-linear-waves-3d test-athena-linear-waves-3d.cpp ${MHD_SOURCES})
deal_ii_setup_target(test-athena-linear-waves-3d)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/mhd/athena-linear-waves)
target_include_directories(test-athena-linear-waves-3d PUBLIC ${CONFIG_HEADER_DIR})
target_include_directories(test-athena-linear-waves-3d PUBLIC
                           ${PROJECT_SOURCE_DIR}/include/sapphirepp/mhd)
target_link_libraries(test-athena-linear-waves-3d UtilsLib)

add_test(NAME athena-linear-waves-3d-sound-wave
         COMMAND ./test-athena-linear-waves-3d parameter.prm 1e-5)
set_property(TEST athena-linear-waves-3d-sound-wave PROPERTY
             LABELS "athena-linear-waves"
                    "3d" "mhd" "sound-wave"
                    "ERK4" "periodic BC"
                    "full-suite-Release" "full-suite-Debug")


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
