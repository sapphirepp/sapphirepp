file(COPY parameter-CN.prm parameter-ERK4.prm parameter-LSERK.prm
          parameter-FE.prm parameter-BE.prm
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


add_executable(test-convergence-study test-convergence-study.cpp ${VFP_SOURCES})
deal_ii_setup_target(test-convergence-study)
file(GLOB CONFIG_HEADER_DIR ${PROJECT_SOURCE_DIR}/examples/convergence-study)
target_include_directories(test-convergence-study PUBLIC ${PROJECT_SOURCE_DIR}/include/sapphirepp/vfp)
target_include_directories(test-convergence-study PUBLIC ${CONFIG_HEADER_DIR})
target_link_libraries(test-convergence-study UtilsLib)

add_test(NAME convergence-study-CN
         COMMAND ./test-convergence-study parameter-CN.prm 5e-6)
set_property(TEST convergence-study-CN PROPERTY
             LABELS "convergence-study"
                    "1d" "spatial_advection" "magnetic"
                    "Crank-Nicholson" "periodic BC"
                    "fast" "commit-Debug" "commit-Release"
                    "full-suite-Release" "full-suite-Debug")

add_test(NAME convergence-study-ERK4
         COMMAND ./test-convergence-study parameter-ERK4.prm 1e-8)
set_property(TEST convergence-study-ERK4 PROPERTY
             LABELS "convergence-study"
                    "1d" "scattering"
                    "ERK4" "periodic BC"
                    "fast" "commit-Debug" "commit-Release"
                    "full-suite-Release" "full-suite-Debug")

# add_test(NAME convergence-study-LSERK
#          COMMAND ./test-convergence-study parameter-LSERK.prm 1e-8)
# set_property(TEST convergence-study-LSERK PROPERTY
#              LABELS "convergence-study"
#                     "1d" "scattering"
#                     "LSERK" "periodic BC"
#                     "full-suite-Debug" "full-suite-Release")

add_test(NAME convergence-study-FE
         COMMAND ./test-convergence-study parameter-FE.prm 2e-3)
set_property(TEST convergence-study-FE PROPERTY
             LABELS "convergence-study"
                    "1d" "scattering"
                    "FE" "periodic BC"
                    "full-suite-Debug" "full-suite-Release")

add_test(NAME convergence-study-BE
         COMMAND ./test-convergence-study parameter-BE.prm 6e-4)
set_property(TEST convergence-study-BE PROPERTY
             LABELS "convergence-study"
                    "1d" "scattering"
                    "BE" "periodic BC"
                    "full-suite-Debug" "full-suite-Release")


set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/results"
)
