file(COPY parameter-1d-CN.prm parameter-1d-ERK4.prm parameter-1d-LSERK.prm
          parameter-1d-FE.prm parameter-1d-BE.prm
          parameter-2d-ERK4.prm
          parameter-3d-ERK4.prm
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


add_executable(test-scattering-only-1d test-scattering-only-1d.cpp ${VFP_SOURCES})
deal_ii_setup_target(test-scattering-only-1d)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/scattering-only)
target_include_directories(test-scattering-only-1d PUBLIC ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
target_include_directories(test-scattering-only-1d PUBLIC ${CONFIG_HEADER_DIR})
target_link_libraries(test-scattering-only-1d UtilsLib)

add_test(NAME scattering-only-1d-CN
         COMMAND ./test-scattering-only-1d parameter-1d-CN.prm test_run 1e-2)
set_property(TEST scattering-only-1d-CN PROPERTY
             LABELS "scattering-only"
                    "1d" "scattering"
                    "CN" "periodic BC"
                    "full-suite-Release" "full-suite-Debug")

add_test(NAME scattering-only-1d-ERK4
         COMMAND ./test-scattering-only-1d parameter-1d-ERK4.prm test_run 1e-4)
set_property(TEST scattering-only-1d-ERK4 PROPERTY
             LABELS "scattering-only"
                    "1d" "scattering"
                    "ERK4" "periodic BC"
                    "fast" "commit-Debug" "commit-Release"
                    "full-suite-Release" "full-suite-Debug")

# add_test(NAME scattering-only-1d-LSERK
#          COMMAND ./test-scattering-only-1d parameter-1d-LSERK.prm test_run 1e-4)
# set_property(TEST scattering-only-1d-LSERK PROPERTY
#              LABELS "scattering-only"
#                     "1d" "scattering"
#                     "LSERK" "periodic BC"
#                     "full-suite-Debug" "full-suite-Release")

add_test(NAME scattering-only-1d-FE
         COMMAND ./test-scattering-only-1d parameter-1d-FE.prm test_run 1e-1)
set_property(TEST scattering-only-1d-FE PROPERTY
             LABELS "scattering-only"
                    "1d" "scattering"
                    "FE" "periodic BC"
                    "full-suite-Debug" "full-suite-Release")

add_test(NAME scattering-only-1d-BE
         COMMAND ./test-scattering-only-1d parameter-1d-BE.prm test_run 1e-1)
set_property(TEST scattering-only-1d-BE PROPERTY
             LABELS "scattering-only"
                    "1d" "scattering"
                    "BE" "periodic BC"
                    "full-suite-Debug" "full-suite-Release")


add_executable(test-scattering-only-2d test-scattering-only-2d.cpp ${VFP_SOURCES})
deal_ii_setup_target(test-scattering-only-2d)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/scattering-only)
target_include_directories(test-scattering-only-2d PUBLIC ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
target_include_directories(test-scattering-only-2d PUBLIC ${CONFIG_HEADER_DIR})
target_link_libraries(test-scattering-only-2d UtilsLib)

add_test(NAME scattering-only-2d
         COMMAND ./test-scattering-only-2d parameter-2d-ERK4.prm 1e-4)
set_property(TEST scattering-only-2d PROPERTY
             LABELS "scattering-only"
                    "2d" "scattering"
                    "ERK4" "periodic BC"
                    "full-suite-Release")


add_executable(test-scattering-only-3d test-scattering-only-3d.cpp ${VFP_SOURCES})
deal_ii_setup_target(test-scattering-only-3d)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/scattering-only)
target_include_directories(test-scattering-only-3d PUBLIC ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
target_include_directories(test-scattering-only-3d PUBLIC ${CONFIG_HEADER_DIR})
target_link_libraries(test-scattering-only-3d UtilsLib)

add_test(NAME scattering-only-3d
         COMMAND ./test-scattering-only-3d parameter-3d-ERK4.prm 1e-4)
set_property(TEST scattering-only-3d PROPERTY
             LABELS "scattering-only"
                    "3d" "scattering"
                    "ERK4" "periodic BC"
                    "commit-Release"
                    "full-suite-Release")


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
