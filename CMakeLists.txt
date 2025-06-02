cmake_minimum_required(VERSION 3.18)
project(astrotiger CXX CUDA)
set(BUILD_OPTIONS_FILE
    "${CMAKE_SOURCE_DIR}/BuildOptions.cmake"
    CACHE FILEPATH
    "Path to a CMake file that sets preprocessor constants")
if(EXISTS "${BUILD_OPTIONS_FILE}")
  message(STATUS "Including build options from ${BUILD_OPTIONS_FILE}")
  include("${BUILD_OPTIONS_FILE}")
endif()
set(DIMENSION_COUNT 3 CACHE STRING "Number of spatial dimensions (must be >= 1)")
if(DIMENSION_COUNT LESS 1)
  message(FATAL_ERROR "DIMENSION_COUNT must be at least 1 (got ${DIMENSION_COUNT})")
endif()
set(MODE_COUNT 3 CACHE STRING "MODE_COUNT of numerical scheme (must be >= 1)")
if(MODE_COUNT LESS 1)
  message(FATAL_ERROR "MODE_COUNT must be at least 1 (got ${MODE_COUNT})")
endif()
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(HPX    REQUIRED)
set(astrotiger_src
  src/main.cpp
  src/Options.cpp
  src/Util.cpp
)
set(astrotiger_header
  include/Array.hpp
  include/Basis.hpp
  include/ContainerArithmetic.hpp
  include/EulerState.hpp
  include/Hdf5.hpp
  include/HyperGrid.hpp
  include/HydroGrid.hpp
  include/Legendre.hpp
  include/Matrix.hpp
  include/MultiIndex.hpp
  include/Octant.hpp
  include/Options.hpp
  include/Quadrature.hpp
  include/Range.hpp
  include/Real.hpp
  include/RungeKutta.hpp
  include/Util.hpp
  include/Zobrist.hpp
)
add_executable(codegen
  src/codegen.cpp
  src/permutation.cpp
  src/Util.cpp
  include/Indent.hpp
  include/Util.hpp
)
target_include_directories(codegen PRIVATE
  ${PROJECT_SOURCE_DIR}/include
)
set(GENERATED_SOURCE_DIR
    ${CMAKE_BINARY_DIR}/generated_source
)
file(MAKE_DIRECTORY ${GENERATED_SOURCE_DIR})
add_custom_target(generate_code
  COMMAND codegen --output-dir ${GENERATED_SOURCE_DIR}
  DEPENDS codegen
  COMMENT "Running codegen → ${GENERATED_SOURCE_DIR}"
)
file(GLOB_RECURSE GENERATED_SOURCES
  CONFIGURE_DEPENDS
  "${GENERATED_SOURCE_DIR}/*.cpp"
)
add_hpx_executable(
  astrotiger
  DEPENDENCIES
    HDF5::HDF5
    siloh5
    stdc++_libbacktrace
  SOURCES
    ${astrotiger_src}
    ${GENERATED_SOURCES}
  HEADERS
    ${astrotiger_header}
)
add_dependencies(astrotiger generate_code)
target_include_directories(astrotiger
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${GENERATED_SOURCE_DIR}
)
add_compile_definitions(
  DIMENSION_COUNT=${DIMENSION_COUNT}
  MODE_COUNT=${MODE_COUNT}
)
set(ALL_FLAGS -std=c++23 -Wall -Wno-comma-subscript -I${PROJECT_SOURCE_DIR}/include)
set(CXX_FLAGS  ${ALL_FLAGS} -march=native)
set(CUDA_FLAGS ${ALL_FLAGS} -Xcompiler -march=native)
set_property(TARGET astrotiger PROPERTY CUDA_ARCHITECTURES native)
set_property(TARGET astrotiger PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_compile_options(astrotiger PUBLIC
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O0 -g -fno-omit-frame-pointer>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O3>>
  $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -g -pg -DNDEBUG>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} -G>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} --use_fast_math -DNDEBUG Xcompiler=-Ofast>>
  $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} -Xcompiler=-g -pg -DNDEBUG>>
)
target_link_options(astrotiger PUBLIC
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -rdynamic>>
)
target_link_libraries(codegen PUBLIC
  stdc++_libbacktrace
)
target_compile_options(codegen PUBLIC
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -O0 -g -fno-omit-frame-pointer>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O2>>
)
target_link_options(codegen PUBLIC
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -rdynamic>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS}>>
  $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:-pg ${CXX_FLAGS}>>
)
