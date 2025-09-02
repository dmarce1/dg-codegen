cmake_minimum_required(VERSION 3.18)

project(astrotiger CXX CUDA)

# ------------------------------------------------------------------------------
# FLAGS
# ------------------------------------------------------------------------------
set(ALL_FLAGS -std=c++23 -Wall -Wno-comma-subscript -I${PROJECT_SOURCE_DIR}/include)
set(CXX_FLAGS  ${ALL_FLAGS} -fext-numeric-literals -march=native)
set(CUDA_FLAGS ${ALL_FLAGS} -Xcompiler -march=native)

# ------------------------------------------------------------------------------
# DEPENDENCIES
# ------------------------------------------------------------------------------
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(HPX  REQUIRED)

# ------------------------------------------------------------------------------
# GOOGLE TEST SETUP
# ------------------------------------------------------------------------------
#include(FetchContent)
#FetchContent_Declare(
#  googletest
#  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
#)
#set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
#FetchContent_MakeAvailable(googletest)

#enable_testing()
#include(GoogleTest)

# ------------------------------------------------------------------------------
# SOURCES
# ------------------------------------------------------------------------------
set(astrotiger_src
  src/main.cpp
  src/Constants.cpp
  src/FermiDirac.cpp
  src/Options.cpp
  src/Util.cpp
  src/Radiation.cpp
  src/SodShockTube.cpp
  src/Symbolic.cpp
)

set(astrotiger_header
   include/AutoDifferentiation.hpp
   include/BiconjugateGradient.hpp
   include/Constants.hpp
   include/ContainerArithmetic.hpp
   include/Definitions.hpp
   include/DoubleReal.hpp
   include/EulerState.hpp
   include/Face.hpp
   include/Hdf5.hpp
   include/HyperSubgrid.hpp
   include/Indent.hpp
   include/Integrate.hpp
   include/Matrix.hpp
   include/MultiIndex.hpp
   include/Octogrid.hpp
   include/Options.hpp
   include/Polynomial.hpp
   include/RadiationState.hpp
   include/Range.hpp
   include/Real.hpp
   include/RungeKutta.hpp
   include/Symbolic.hpp
   include/TriIndex.hpp
   include/Util.hpp
)


# ------------------------------------------------------------------------------
# CODEGEN
# ------------------------------------------------------------------------------
add_executable(codegen
  	src/codegen.cpp
  	src/stateCodegen.cpp
  	src/Util.cpp
  	include/Definitions.hpp
  	include/Indent.hpp
  	include/stateCodegen.hpp
	include/Symbolic.hpp
 	include/Util.hpp
)

target_include_directories(codegen PRIVATE
  ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(codegen PUBLIC
    gmp
	stdc++_libbacktrace
    quadmath    
    symengine
)

target_compile_options(codegen PUBLIC
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -O0 -g -fno-omit-frame-pointer>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -O2>>
  $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -g -pg -march=native -O2 -DNDEBUG>>
)

target_link_options(codegen PUBLIC
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native -rdynamic>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -march=native>>
  $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -pg -march=native>>
)

# ------------------------------------------------------------------------------
# GENERATED CODE
# ------------------------------------------------------------------------------
set(GENERATED_SOURCE_DIR ${CMAKE_BINARY_DIR}/generated_source)
file(MAKE_DIRECTORY ${GENERATED_SOURCE_DIR})
set(GENERATED_DG_TRANSFORMS "${GENERATED_SOURCE_DIR}/dgTransforms.hpp")

add_custom_command(
  OUTPUT  "${GENERATED_DG_TRANSFORMS}"
          "${GENERATED_SOURCE_DIR}/.codegen.stamp"
  COMMAND "${CMAKE_COMMAND}" -E make_directory "${GENERATED_SOURCE_DIR}"
  COMMAND codegen --output-dir "${GENERATED_SOURCE_DIR}"
  COMMAND "${CMAKE_COMMAND}" -E compare_files
          "${GENERATED_DG_TRANSFORMS}" "${GENERATED_DG_TRANSFORMS}"
  COMMAND "${CMAKE_COMMAND}" -E touch
          "${GENERATED_SOURCE_DIR}/.codegen.stamp"
  DEPENDS codegen
  COMMENT "Running codegen → ${GENERATED_SOURCE_DIR}"
  VERBATIM
)

add_custom_target(generate_code
  DEPENDS "${GENERATED_DG_TRANSFORMS}"
)

file(GLOB_RECURSE GENERATED_SOURCES
  CONFIGURE_DEPENDS
  "${GENERATED_SOURCE_DIR}/*.cpp"
)

# ------------------------------------------------------------------------------
# ASTROTIGER EXECUTABLE
# ------------------------------------------------------------------------------
add_hpx_executable(
  astrotiger
  DEPENDENCIES
    HDF5::HDF5
    stdc++_libbacktrace
    quadmath
    mpfr
    gmp  
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

set_property(TARGET astrotiger PROPERTY CUDA_ARCHITECTURES native)
set_property(TARGET astrotiger PROPERTY CUDA_SEPARABLE_COMPILATION ON)

target_compile_options(astrotiger PUBLIC
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O3 -g -fopt-info-vec-optimized >>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -O0 -g -fno-omit-frame-pointer>>
  $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -g -pg -O2 -DNDEBUG>>
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} -G>>
  $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} --use_fast_math -DNDEBUG Xcompiler=-Ofast>>
  $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS} -Xcompiler=-g -pg -DNDEBUG>>
)

target_link_options(astrotiger PUBLIC
  $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS} -rdynamic>>
)

