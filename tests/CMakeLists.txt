include(GoogleTest)

set(FG_UNIT_TESTS
  test_index_map
  test_coefficients
  test_toeplitz
  test_smatrix
  test_butcher
  test_propagation_matrix
  test_incidence
)
foreach(t ${FG_UNIT_TESTS})
  add_executable(${t} unit/${t}.cpp)
  target_link_libraries(${t} PRIVATE fg_core GTest::gtest_main)
  gtest_discover_tests(${t} PROPERTIES LABELS "unit")
endforeach()

# Validation cases (Plan section "Validierung") are added here from milestone M2 on, e.g.
#   add_test(NAME V1_fresnel COMMAND fgsolve ${PROJECT_SOURCE_DIR}/examples/v1_fresnel.json)
#   set_tests_properties(V1_fresnel PROPERTIES LABELS "validation")
