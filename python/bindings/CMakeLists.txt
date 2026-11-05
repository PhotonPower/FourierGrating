# Python extension module <FG_PYTHON_PACKAGE>._core (Plan section "M12: Python-Anbindung").
pybind11_add_module(_core MODULE module.cpp)
target_link_libraries(_core PRIVATE fg_core nlohmann_json::nlohmann_json fg_warnings)
target_compile_definitions(_core PRIVATE FG_PYTHON_PACKAGE="${FG_PYTHON_PACKAGE}")
# only this component goes into the wheel (Eigen, nlohmann_json etc. install headers otherwise)
install(TARGETS _core LIBRARY DESTINATION ${FG_PYTHON_PACKAGE} COMPONENT python)
