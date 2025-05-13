set_directory_properties(PROPERTIES LABELS "mhd")

# add_subdirectory(athena-linear-waves)
add_subdirectory(linear-waves)
add_subdirectory(shock-tube)


add_executable(test-mhd-equations test-mhd-equations.cpp)
deal_ii_setup_target(test-mhd-equations)
target_link_libraries(test-mhd-equations UtilsLib MHDLib)

add_test(NAME test-mhd-equations COMMAND
  ./test-mhd-equations
)
set_property(TEST test-mhd-equations PROPERTY LABELS
  "mhd-equations"
  "unit-test"
  "fast" "commit-Debug" "commit-Release"
  "full-suite-Release" "full-suite-Debug"
)
