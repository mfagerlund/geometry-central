cmake_minimum_required(VERSION 3.8.0)

project(geodesic-comparison)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

### Configure the compiler
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++14 -Wall -Wextra -g3")
  SET(CMAKE_CXX_FLAGS_DEBUG          "")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized")
  endif()

  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  SET(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -g3")
  SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions (-DNOMINMAX)
else()
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()

# Build the benchmarks
add_executable(geodesic-comparison geodesic_comparison.cpp)
target_link_libraries(geodesic-comparison geometry-central)

add_executable(gfr-strategy-benchmark gfr_strategy_benchmark.cpp)
target_link_libraries(gfr-strategy-benchmark geometry-central)

add_executable(multi-mesh-benchmark multi_mesh_benchmark.cpp)
target_link_libraries(multi-mesh-benchmark geometry-central)

add_executable(debug-path debug_path.cpp)
target_link_libraries(debug-path geometry-central)

add_executable(initial-distance-test initial_distance_test.cpp)
target_link_libraries(initial-distance-test geometry-central)

add_executable(diagnose-facestrip diagnose_facestrip.cpp)
target_link_libraries(diagnose-facestrip geometry-central)

add_executable(fan-direction-test fan_direction_test.cpp)
target_link_libraries(fan-direction-test geometry-central)

add_executable(cold-cache-benchmark cold_cache_benchmark.cpp)
target_link_libraries(cold-cache-benchmark geometry-central)

add_executable(epsilon-test epsilon_test.cpp)
target_link_libraries(epsilon-test geometry-central)

add_executable(epsilon-vs-flipout-test epsilon_vs_flipout_test.cpp)
target_link_libraries(epsilon-vs-flipout-test geometry-central)

# Add geometry central as a subproject
add_subdirectory(../ "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
