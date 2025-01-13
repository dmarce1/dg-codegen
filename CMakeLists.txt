cmake_minimum_required(VERSION 3.18)

project(astrotiger CXX CUDA)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(HPX REQUIRED)

set(AT_NDIM "3" CACHE STRING "Number of spatial dimensions")
set(AT_LEGENDRE_MAX_ORDER "32" CACHE STRING "Sets maximum order for Legendre polynomials")
set(AT_CHECK_REALS "1" CACHE STRING "Whether to check Real number operations for errors.")

set(ALL_FLAGS  -Wall -I${PROJECT_SOURCE_DIR}/include -DNDIM=${AT_NDIM} -DCHECK_REALS=${AT_CHECK_REALS} -DLEGENDRE_MAX_ORDER=${AT_LEGENDRE_MAX_ORDER})
set(CXX_FLAGS  ${ALL_FLAGS} -march=native )
set(CUDA_FLAGS ${ALL_FLAGS} -Xcompiler -march=native )

set(astrotiger_src
	src/main.cpp
	src/Options.cpp
	src/Sod.cpp
) 

set(astrotiger_header
	include/Complex.hpp
	include/Sod.hpp
	include/Integrate.hpp
	include/Matrix.hpp
	include/Numbers.hpp
	include/LegendreP.hpp
	include/Options.hpp
	include/Polynomial.hpp
	include/Real.hpp
	include/Root.hpp
	include/Quadrature.hpp
	include/Utilities.hpp
	include/Vector.hpp
) 

add_hpx_executable(
  astrotiger
  DEPENDENCIES
  stdc++_libbacktrace
  SOURCES
    ${astrotiger_src}
  HEADERS
    ${astrotiger_header}
)
set_property(TARGET astrotiger PROPERTY CUDA_ARCHITECTURES native)
set_property(TARGET astrotiger PROPERTY CUDA_SEPARABLE_COMPILATION ON)

target_link_options(
	astrotiger
	PUBLIC
		$<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:
			${CXX_FLAGS}
   	    >>
		$<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:
			-pg ${CXX_FLAGS}
   	    >>
)

target_compile_options(
    astrotiger
    PUBLIC
        $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CXX>:${CXX_FLAGS}>>
        $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CUDA>:-G ${CUDA_FLAGS}>>
        
        $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:-Ofast ${CXX_FLAGS}>>
        $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CUDA>:
            --use_fast_math -Xcompiler=-Ofast ${CUDA_FLAGS}>>
        
        $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CXX>:-O2 -pg ${CXX_FLAGS}>>
        $<$<CONFIG:Profile>:$<$<COMPILE_LANGUAGE:CUDA>:
            -Xcompiler=-O2 -pg ${CUDA_FLAGS}>>
)

