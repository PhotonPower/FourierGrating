# Dependency resolution: prefer system packages (AlmaLinux 10: dnf), fall back to FetchContent.
include(FetchContent)

if(MSVC)
  # MSVC (M12): never pick up MinGW packages from an MSYS2 installation on the PATH (their headers
  # do not compile with cl), nor copies registered by other projects in the CMake package registry.
  list(APPEND CMAKE_IGNORE_PREFIX_PATH C:/msys64/ucrt64 C:/msys64/mingw64 C:/msys64/clang64 C:/msys64/usr)
  set(CMAKE_FIND_USE_PACKAGE_REGISTRY OFF)
endif()

# --- Eigen3 (header-only linear algebra, CPU reference backend) ---
find_package(Eigen3 3.4 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  message(STATUS "Eigen3 not found on the system - fetching 3.4.0 (headers only)")
  # SOURCE_SUBDIR without CMakeLists.txt: Eigen's own CMake (BLAS/Fortran probes, which fail with a
  # foreign gfortran on the PATH, e.g. GitHub's Windows runners) is not run; Eigen is header-only.
  FetchContent_Declare(eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    SOURCE_SUBDIR fg-headers-only)
  FetchContent_MakeAvailable(eigen)
  add_library(fg_eigen INTERFACE)
  target_include_directories(fg_eigen SYSTEM INTERFACE ${eigen_SOURCE_DIR})
  add_library(Eigen3::Eigen ALIAS fg_eigen)
endif()

# --- nlohmann_json (configuration files) ---
find_package(nlohmann_json 3.10 QUIET)
if(NOT TARGET nlohmann_json::nlohmann_json)
  message(STATUS "nlohmann_json not found on the system - fetching v3.11.3")
  FetchContent_Declare(nlohmann_json
    URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  FetchContent_MakeAvailable(nlohmann_json)
endif()

# --- OpenMP (optional) ---
if(FG_ENABLE_OPENMP)
  find_package(OpenMP QUIET)
endif()

# --- HDF5 (optional) ---
if(FG_ENABLE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS CXX)
endif()

# --- pybind11 (Python module, M12; scikit-build-core provides it in the build environment) ---
if(FG_BUILD_PYTHON)
  find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 2.12 CONFIG QUIET)
  if(NOT pybind11_FOUND)
    message(STATUS "pybind11 not found - fetching v2.13.6")
    FetchContent_Declare(pybind11
      URL https://github.com/pybind/pybind11/archive/refs/tags/v2.13.6.tar.gz
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    FetchContent_MakeAvailable(pybind11)
  endif()
endif()

# --- GoogleTest ---
if(FG_BUILD_TESTS)
  find_package(GTest QUIET)
  if(NOT TARGET GTest::gtest_main)
    message(STATUS "GoogleTest not found on the system - fetching v1.15.2")
    FetchContent_Declare(googletest
      URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
  endif()
endif()
