configure_file(${PROJECT_SOURCE_DIR}/include/fg/version.hpp.in
               ${PROJECT_BINARY_DIR}/generated/fg/version.hpp @ONLY)

add_library(fg_core
  core/index_map.cpp
  core/quadrature.cpp
  fourier/coefficients.cpp
  fourier/toeplitz.cpp
  fourier/sampling.cpp
  fourier/contour.cpp
  geometry/stack.cpp
  geometry/contour.cpp
  geometry/normal_field.cpp
  geometry/profile.cpp
  physics/propagation_matrix.cpp
  physics/modes.cpp
  physics/incidence.cpp
  physics/qmatrix.cpp
  physics/efficiency.cpp
  physics/anisotropic.cpp
  smatrix/smatrix.cpp
  solvers/butcher.cpp
  solvers/solver.cpp
  solvers/layer_model.cpp
  solvers/rcwa.cpp
  solvers/rk.cpp
  solvers/simulation.cpp
  field/nearfield.cpp
  linalg/cpu_backend.cpp
  io/config.cpp
)
add_library(FourierGrating::core ALIAS fg_core)

target_include_directories(fg_core PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/generated>)
target_link_libraries(fg_core PUBLIC Eigen3::Eigen PRIVATE nlohmann_json::nlohmann_json fg_warnings)
target_compile_features(fg_core PUBLIC cxx_std_20)

if(FG_ENABLE_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(fg_core PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(fg_core PUBLIC FG_HAVE_OPENMP)
endif()

if(FG_ENABLE_HDF5)
  target_link_libraries(fg_core PUBLIC HDF5::HDF5)
  target_compile_definitions(fg_core PUBLIC FG_HAVE_HDF5)
endif()

if(FG_ENABLE_CUDA)
  target_sources(fg_core PRIVATE cuda/cuda_backend.cpp)
  target_link_libraries(fg_core PUBLIC CUDA::cudart CUDA::cublas CUDA::cusolver)
  target_compile_definitions(fg_core PUBLIC FG_HAVE_CUDA)
  if(MINGW)
    # See cmake/mingw_cuda_security_stub.cpp: the import libs above need MSVC CRT symbols that
    # MinGW lacks. Deliberately NOT linked into fg_core: fg_core is a static library, and CMake
    # merges an OBJECT library linked to a static library into its archive -- which makes the
    # stub an archive member again, only pulled in if ld already needs it by the time it scans
    # fg_core.a, i.e. before it reaches cudart.lib later on the link line. It must instead reach
    # each final executable's own link line as a plain object (see top-level CMakeLists.txt).
    add_library(fg_cuda_mingw_abi_shim OBJECT ${PROJECT_SOURCE_DIR}/cmake/mingw_cuda_security_stub.cpp)
  endif()
endif()
