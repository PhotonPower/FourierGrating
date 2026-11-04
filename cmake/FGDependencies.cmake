# Dependency resolution: prefer system packages (AlmaLinux 10: dnf), fall back to FetchContent.
include(FetchContent)

# --- Eigen3 (header-only linear algebra, CPU reference backend) ---
find_package(Eigen3 3.4 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  message(STATUS "Eigen3 not found on the system - fetching 3.4.0")
  FetchContent_Declare(eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(eigen)
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
