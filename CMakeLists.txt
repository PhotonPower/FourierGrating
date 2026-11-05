cmake_minimum_required(VERSION 3.24)
project(FourierGrating VERSION 0.14.0 DESCRIPTION "Fourier modal method: differential method (implicit RK) and RCWA" LANGUAGES CXX)

option(FG_ENABLE_CUDA   "Build the CUDA backend (RTX 6000 Ada, sm_89)" OFF)
option(FG_ENABLE_OPENMP "Use OpenMP for parameter sweeps"              ON)
option(FG_ENABLE_HDF5   "Enable HDF5 output for near fields"           OFF)
option(FG_BUILD_TESTS   "Build unit and validation tests"              ON)
option(FG_BUILD_APPS    "Build command line applications"              ON)
option(FG_WARNINGS_AS_ERRORS "Treat compiler warnings as errors"       OFF)
option(FG_BUILD_PYTHON  "Build the Python module (scikit-build sets ON)" OFF)
option(FG_ENABLE_NURBFORGE "CAD slicer: STEP files via Nurbforge (layer type cad, M14)" OFF)
set(FG_NURBFORGE_TAG "v0.6.0" CACHE STRING "Nurbforge tag fetched with FG_ENABLE_NURBFORGE")
set(FG_PYTHON_PACKAGE "fouriergrating" CACHE STRING
    "Python package receiving the _core module (fouriergrating | fouriergrating_cuda12)")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FGCompilerFlags)
include(FGDependencies)

if(FG_ENABLE_CUDA)
  # M10: host API of cuBLAS/cuSOLVER only (no device kernels yet) -> no CUDA language/nvcc needed.
  # CMAKE_CUDA_ARCHITECTURES (89, RTX 6000 Ada) becomes relevant once own kernels are added.
  if(MINGW AND DEFINED ENV{CUDA_PATH})
    # FindCUDAToolkit only picks the "lib/x64" search suffix for an MSVC-style architecture id
    # (AMD64/x64); GCC reports "x86_64", so it never finds cudart.lib without this hint on
    # Windows/MinGW (see docs/build_windows.md).
    list(APPEND CMAKE_LIBRARY_PATH "$ENV{CUDA_PATH}/lib/x64")
  endif()
  find_package(CUDAToolkit 12.6 REQUIRED)   # cusolverDnXgeev requires CUDA >= 12.6 Update 2
endif()

if(FG_BUILD_PYTHON)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # fg_core is linked into a shared Python module
endif()

add_subdirectory(src)

if(FG_ENABLE_CUDA AND MINGW)
  # fg_cuda_mingw_abi_shim (src/CMakeLists.txt) must reach every executable's own link line as a
  # plain object, not fg_core's archive -- see the comment there. link_libraries() applies to all
  # targets defined from here on (apps, tests), which are executables, not static libraries.
  link_libraries(fg_cuda_mingw_abi_shim)
endif()

if(FG_BUILD_APPS)
  add_subdirectory(apps)
endif()
if(FG_BUILD_PYTHON)
  add_subdirectory(python/bindings)
endif()
if(FG_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

message(STATUS "FourierGrating ${PROJECT_VERSION}: CUDA=${FG_ENABLE_CUDA} OpenMP=${FG_ENABLE_OPENMP} HDF5=${FG_ENABLE_HDF5} Tests=${FG_BUILD_TESTS} Python=${FG_BUILD_PYTHON} Nurbforge=${FG_ENABLE_NURBFORGE}")
