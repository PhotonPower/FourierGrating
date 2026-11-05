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
  test_profile
  test_fff_pairing
  test_anisotropic
  test_normal_field2d
  test_backend
  test_scan
  test_polarization_split
  test_contour
  test_normal_field_samples
)
foreach(t ${FG_UNIT_TESTS})
  add_executable(${t} unit/${t}.cpp)
  target_link_libraries(${t} PRIVATE fg_core GTest::gtest_main)
  gtest_discover_tests(${t} PROPERTIES LABELS "unit")
endforeach()
add_executable(test_simulation unit/test_simulation.cpp)
target_link_libraries(test_simulation PRIVATE fg_core GTest::gtest_main nlohmann_json::nlohmann_json)
gtest_discover_tests(test_simulation PROPERTIES LABELS "unit")

# Validation cases (Plan section "Validierung") are added here from milestone M2 on, e.g.
#   add_test(NAME V1_fresnel COMMAND fgsolve ${PROJECT_SOURCE_DIR}/examples/v1_fresnel.json)
#   set_tests_properties(V1_fresnel PROPERTIES LABELS "validation")

# ---- validation cases (label "validation") ----
add_executable(test_validation_v1_v3 validation/test_v1_v3.cpp)
target_link_libraries(test_validation_v1_v3 PRIVATE fg_core GTest::gtest_main)
gtest_discover_tests(test_validation_v1_v3 PROPERTIES LABELS "validation")
add_executable(test_validation_v4_v6 validation/test_v4_v6.cpp)
target_link_libraries(test_validation_v4_v6 PRIVATE fg_core GTest::gtest_main)
gtest_discover_tests(test_validation_v4_v6 PROPERTIES LABELS "validation" TIMEOUT 900 DISCOVERY_TIMEOUT 60)
add_executable(test_validation_v7_v8 validation/test_v7_v8.cpp)
target_link_libraries(test_validation_v7_v8 PRIVATE fg_core GTest::gtest_main)
gtest_discover_tests(test_validation_v7_v8 PROPERTIES LABELS "validation" TIMEOUT 900 DISCOVERY_TIMEOUT 60)
add_executable(test_validation_v12 validation/test_v12_nearfield.cpp)
target_link_libraries(test_validation_v12 PRIVATE fg_core GTest::gtest_main nlohmann_json::nlohmann_json)
gtest_discover_tests(test_validation_v12 PROPERTIES LABELS "validation" TIMEOUT 900 DISCOVERY_TIMEOUT 60)
add_executable(test_validation_v13_v14 validation/test_v13_gold.cpp)
target_link_libraries(test_validation_v13_v14 PRIVATE fg_core GTest::gtest_main)
gtest_discover_tests(test_validation_v13_v14 PROPERTIES LABELS "validation" TIMEOUT 900 DISCOVERY_TIMEOUT 60)

# ---- command line smoke tests on example configurations ----
foreach(ex v1_fresnel_interface v3_lamellar_dielectric v3_lamellar_rk_scan v5_glass_deep_te echelette_au_scan v7_cobalt_anisotropic nearfield_lamellar_1d nearfield_pillars_2d sweep_lamellar_depth)
  add_test(NAME example_${ex}
           COMMAND fgsolve ${PROJECT_SOURCE_DIR}/examples/${ex}.json
           WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  set_tests_properties(example_${ex} PROPERTIES LABELS "validation")
endforeach()
add_executable(test_validation_v9_v10 validation/test_v9_v10.cpp)
target_link_libraries(test_validation_v9_v10 PRIVATE fg_core GTest::gtest_main)
target_compile_definitions(test_validation_v9_v10 PRIVATE FG_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
gtest_discover_tests(test_validation_v9_v10 PROPERTIES LABELS "validation" TIMEOUT 900 DISCOVERY_TIMEOUT 60)
add_executable(test_validation_v8_2d validation/test_v8_2d.cpp)
target_link_libraries(test_validation_v8_2d PRIVATE fg_core GTest::gtest_main)
gtest_discover_tests(test_validation_v8_2d PROPERTIES LABELS "validation" TIMEOUT 900 DISCOVERY_TIMEOUT 60)
add_executable(test_validation_v15 validation/test_v15_contour.cpp)
target_link_libraries(test_validation_v15 PRIVATE fg_core GTest::gtest_main)
gtest_discover_tests(test_validation_v15 PROPERTIES LABELS "validation" TIMEOUT 900 DISCOVERY_TIMEOUT 60)
add_executable(test_validation_v16a validation/test_v16a_sectioned.cpp)
target_link_libraries(test_validation_v16a PRIVATE fg_core GTest::gtest_main)
gtest_discover_tests(test_validation_v16a PROPERTIES LABELS "validation" TIMEOUT 1800 DISCOVERY_TIMEOUT 60)
add_executable(test_validation_v12_2d validation/test_v12_2d.cpp)
target_link_libraries(test_validation_v12_2d PRIVATE fg_core GTest::gtest_main nlohmann_json::nlohmann_json)
gtest_discover_tests(test_validation_v12_2d PROPERTIES LABELS "validation" TIMEOUT 900 DISCOVERY_TIMEOUT 60)

# ---- regression protection (version 2): all examples against golden files (label "regression") ----
add_executable(test_regression regression/test_regression.cpp)
target_link_libraries(test_regression PRIVATE fg_core GTest::gtest_main)
target_compile_definitions(test_regression PRIVATE FG_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
# the test list comes from examples.txt: rebuild (and rediscover) when it changes
set_source_files_properties(regression/test_regression.cpp PROPERTIES
  OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/regression/examples.txt)
gtest_discover_tests(test_regression PROPERTIES LABELS "regression" TIMEOUT 1800 DISCOVERY_TIMEOUT 60)
