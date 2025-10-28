cmake_minimum_required(VERSION 3.18)
project(astrotiger C CXX)
set(CMAKE_CXX_STANDARD 23)
add_executable(codegen
  	src/codegen.cpp
   	src/Util.cpp
  	src/Indent.hpp
 	src/Util.hpp
)
target_link_libraries(codegen PUBLIC
    gmp
	stdc++_libbacktrace
    quadmath    
    symengine
)


