cmake_minimum_required(VERSION 3.18)
project(astrotiger CXX CUDA)
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(HPX REQUIRED)
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
add_hpx_executable(
  astrotiger
  DEPENDENCIES
  	siloh5
	stdc++_libbacktrace
  SOURCES
    ${astrotiger_src}
  HEADERS
    ${astrotiger_header}
)
set(ALL_FLAGS -Wall -Wno-comma-subscript -I${PROJECT_SOURCE_DIR}/include )
set(CXX_FLAGS  ${ALL_FLAGS} -march=native )
set(CUDA_FLAGS ${ALL_FLAGS} -Xcompiler -march=native )
set_property(TARGET astrotiger PROPERTY CUDA_ARCHITECTURES native)
set_property(TARGET astrotiger PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_compile_options(
    astrotiger
    PUBLIC
    $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>: ${CXX_FLAGS} -O0 -g -fno-omit-frame-pointer >>
    $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>: ${CXX_FLAGS} -O3 >>
    $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>: ${CXX_FLAGS} -g -pg -DNDEBUG>>
    $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CUDA>: ${CUDA_FLAGS} -G >>
    $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CUDA>: ${CUDA_FLAGS} --use_fast_math -DNDEBUG Xcompiler=-Ofast >>
    $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CUDA>: ${CUDA_FLAGS} -Xcompiler=-g -pg -DNDEBUG>>
)
target_link_options(
	astrotiger
  PUBLIC
	$<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>: ${CXX_FLAGS} -rdynamic>>
	$<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>: ${CXX_FLAGS}>>
   	$<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>: -pg ${CXX_FLAGS}>> 
)
target_link_libraries(
	astrotiger 
  PUBLIC 
	HDF5::HDF5
)
target_link_libraries(
	astrotiger 
  PUBLIC
  	$<$<CONFIG:Debug>:stdc++_libbacktrace>
)



