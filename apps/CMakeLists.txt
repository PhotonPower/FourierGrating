add_executable(fgsolve fgsolve/main.cpp)
target_link_libraries(fgsolve PRIVATE fg_core fg_warnings)
add_executable(fgbench fgbench/main.cpp)
target_link_libraries(fgbench PRIVATE fg_core fg_warnings)
