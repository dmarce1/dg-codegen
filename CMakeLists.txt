cmake_minimum_required(VERSION 3.18)

project(astrotiger CXX CUDA)

set(
   BUILD_OPTIONS_FILE
      "${CMAKE_SOURCE_DIR}/BuildOptions.cmake"
   CACHE FILEPATH
      "Path to a CMake file that sets preprocessor constants"
)

set(ALL_FLAGS -std=c++23 -Wall -Wno-comma-subscript -I${PROJECT_SOURCE_DIR}/include)
set(CXX_FLAGS  ${ALL_FLAGS} -march=native)
set(CUDA_FLAGS ${ALL_FLAGS} -Xcompiler -march=native)
set(DIMENSION_COUNT 3 CACHE STRING "Number of spatial dimensions (must be >= 1)")
set(SIMD_WIDTH 4 CACHE STRING "Width of the native simd type (must be >= 1)")
set(MODE_COUNT 3 CACHE STRING "Number of modes (must be >= 1)")

if(EXISTS "${BUILD_OPTIONS_FILE}")
    message(STATUS "Including build options from ${BUILD_OPTIONS_FILE}")
    include("${BUILD_OPTIONS_FILE}")
endif()
if(DIMENSION_COUNT LESS 1)
  message(FATAL_ERROR "DIMENSION_COUNT must be at least 1 (got ${DIMENSION_COUNT})")
endif()
if(SIMD_WIDTH LESS 1)
  message(FATAL_ERROR "SIMD_WIDTH must be at least 1 (got ${SIMD_WIDTH})")
endif()
if(MODE_COUNT LESS 1)
  message(FATAL_ERROR "MODE_COUNT must be at least 1 (got ${MODE_COUNT})")
endif()

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(HPX  REQUIRED)

set(astrotiger_src
  src/main.cpp
  src/Options.cpp
  src/Util.cpp
)

set(astrotiger_header
  include/Array.hpp
  include/Basis.hpp
  include/ContainerArithmetic.hpp
  include/Definitions.hpp
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

# ---------------------------------------------------------------------
# 1) Build “codegen” exactly as before.  It only rebuilds when its own
#    C++ sources change, because it’s a normal add_executable(...)
# ---------------------------------------------------------------------
add_executable(codegen
  src/codegen.cpp
  src/permutation.cpp
  src/Util.cpp
  include/Definitions.hpp
  include/Indent.hpp
  include/Util.hpp
)

target_include_directories(codegen PRIVATE
  ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(codegen PUBLIC
   stdc++_libbacktrace
)

target_compile_options(codegen PUBLIC
   $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -O0 -g -fno-omit-frame-pointer -fsanitize=address>>
   $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -O2>>
)

target_link_options(codegen PUBLIC
   $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -rdynamic -fsanitize=address>>
   $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native>>
   $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:-pg ${CXX_FLAGS} -march=native>>
)

# ---------------------------------------------------------------------
# 2) Define where generated files will live, and make sure the directory exists
# ---------------------------------------------------------------------
set(GENERATED_SOURCE_DIR
    ${CMAKE_BINARY_DIR}/generated_source
)

# The directory must exist before we run codegen
file(MAKE_DIRECTORY ${GENERATED_SOURCE_DIR})

# ---------------------------------------------------------------------
# 3) Instead of add_custom_target(generate_code ... ALWAYS_RUN),
#    we now use add_custom_command() that produces a “.stamp” file.
#    - The OUTPUT is ${GENERATED_STAMP}
#    - We DEPEND on the codegen executable.
#    - If codegen is rebuilt (timestamp changes) or if ${GENERATED_STAMP} is missing,
#      CMake will rerun the command; otherwise it is up-to-date.
# ---------------------------------------------------------------------
set(GENERATED_STAMP "${GENERATED_SOURCE_DIR}/.codegen.stamp")

add_custom_command(
  OUTPUT ${GENERATED_STAMP}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${GENERATED_SOURCE_DIR}
  COMMAND codegen --output-dir ${GENERATED_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} -E touch ${GENERATED_STAMP}
  DEPENDS codegen
  COMMENT "Running codegen → ${GENERATED_SOURCE_DIR}"
  VERBATIM
)

# Now wrap that single OUTPUT-stamp into a custom target
add_custom_target(generate_code
  DEPENDS ${GENERATED_STAMP}
)

# ---------------------------------------------------------------------
# 4) Glob all generated .cpp files (CMake will re-evaluate this list
#    whenever files under generated_source/ appear or disappear).
#    The CONFIGURE_DEPENDS flag makes CMake rescan the directory
#    at build time, but note that CMake only discovers new files after
#    a previous generation has run. 
# ---------------------------------------------------------------------
file(GLOB_RECURSE GENERATED_SOURCES
  CONFIGURE_DEPENDS
  "${GENERATED_SOURCE_DIR}/*.cpp"
)

# ---------------------------------------------------------------------
# 5) Build the actual “astrotiger” executable.
#    We include both the hand-written src and any generated sources.
#    We also add a dependency so that “generate_code” (and thus codegen)
#    must run (if out of date) before building astrotiger.
# ---------------------------------------------------------------------
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

# Include both your normal headers and the generated directory
target_include_directories(astrotiger
   PUBLIC
      ${PROJECT_SOURCE_DIR}/include
      ${GENERATED_SOURCE_DIR}
)

add_compile_definitions(
   DIMENSION_COUNT=${DIMENSION_COUNT}
   MODE_COUNT=${MODE_COUNT}
   SIMD_WIDTH=${SIMD_WIDTH}
)

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

