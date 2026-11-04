cmake_minimum_required(VERSION 3.24)
project(FourierGrating VERSION 0.8.0 DESCRIPTION "Fourier modal method: differential method (implicit RK) and RCWA" LANGUAGES CXX)

option(FG_ENABLE_CUDA   "Build the CUDA backend (RTX 6000 Ada, sm_89)" OFF)
option(FG_ENABLE_OPENMP "Use OpenMP for parameter sweeps"              ON)
option(FG_ENABLE_HDF5   "Enable HDF5 output for near fields"           OFF)
option(FG_BUILD_TESTS   "Build unit and validation tests"              ON)
option(FG_BUILD_APPS    "Build command line applications"              ON)
option(FG_WARNINGS_AS_ERRORS "Treat compiler warnings as errors"       OFF)

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
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 89)   # RTX 6000 Ada
  endif()
  enable_language(CUDA)
  find_package(CUDAToolkit 12.6 REQUIRED)   # cusolverDnXgeev requires CUDA >= 12.6 Update 2
endif()

add_subdirectory(src)
if(FG_BUILD_APPS)
  add_subdirectory(apps)
endif()
if(FG_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

message(STATUS "FourierGrating ${PROJECT_VERSION}: CUDA=${FG_ENABLE_CUDA} OpenMP=${FG_ENABLE_OPENMP} HDF5=${FG_ENABLE_HDF5} Tests=${FG_BUILD_TESTS}")
