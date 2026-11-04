include(GoogleTest)

set(FG_UNIT_TESTS
  test_index_map
  test_coefficients
  test_toeplitz
  test_smatrix
  test_butcher
  test_propagation_matrix
  test_incidence
  test_qmatrix
  test_rcwa_rk
)
foreach(t ${FG_UNIT_TESTS})
  add_executable(${t} unit/${t}.cpp)
  target_link_libraries(${t} PRIVATE fg_core GTest::gtest_main)
  gtest_discover_tests(${t} PROPERTIES LABELS "unit")
endforeach()

# Validation cases (Plan section "Validierung") are added here from milestone M2 on, e.g.
#   add_test(NAME V1_fresnel COMMAND fgsolve ${PROJECT_SOURCE_DIR}/examples/v1_fresnel.json)
#   set_tests_properties(V1_fresnel PROPERTIES LABELS "validation")

# ---- validation cases (label "validation") ----
add_executable(test_validation_v1_v3 validation/test_v1_v3.cpp)
target_link_libraries(test_validation_v1_v3 PRIVATE fg_core GTest::gtest_main)
gtest_discover_tests(test_validation_v1_v3 PROPERTIES LABELS "validation")

# ---- command line smoke tests on example configurations ----
foreach(ex v1_fresnel_interface v3_lamellar_dielectric v3_lamellar_rk_scan)
  add_test(NAME example_${ex}
           COMMAND fgsolve ${PROJECT_SOURCE_DIR}/examples/${ex}.json
           WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  set_tests_properties(example_${ex} PROPERTIES LABELS "validation")
endforeach()
