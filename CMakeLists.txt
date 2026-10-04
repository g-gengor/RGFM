cmake_minimum_required(VERSION 3.16)
project(RGFM LANGUAGES CXX)

# Default to an optimized build. Release uses -O2 (as the original build
# script did) instead of CMake's default -O3.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()
string(REPLACE "-O3" "-O2" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

# Strict ISO C++14 (no GNU extensions): the sources use the std::complex
# literal `1i`, which GNU dialects would parse as a GCC builtin instead.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Dependencies (all header-only apart from the system thread library).
# Eigen: tested with 3.4 and 5.0, so no version is pinned here.
find_package(Eigen3 REQUIRED NO_MODULE)
# Boost: only Boost.Math headers are used (spherical harmonics and the
# barycentric rational interpolator). Prefer BoostConfig.cmake (Boost >= 1.70);
# fall back to CMake's FindBoost for plain header trees (-DBoost_ROOT=...).
find_package(Boost CONFIG QUIET)
if(NOT Boost_FOUND)
  if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)  # keep FindBoost available on CMake >= 3.30
  endif()
  find_package(Boost REQUIRED MODULE)
endif()
find_package(Threads REQUIRED)

add_executable(run_RGFM
  src/run_RGFM.cpp
  src/RGFM_parser.cpp
  src/Greens_PWexpansion.cpp
)
target_include_directories(run_RGFM PRIVATE src)
target_link_libraries(run_RGFM PRIVATE Eigen3::Eigen Boost::headers Threads::Threads)

# Put the executable in <repo>/bin for every generator ($<1:...> stops
# multi-config generators from appending a per-config subdirectory).
set_target_properties(run_RGFM PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "$<1:${PROJECT_SOURCE_DIR}/bin>"
)
