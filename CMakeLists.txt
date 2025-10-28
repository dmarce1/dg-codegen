cmake_minimum_required(VERSION 3.18)
project(astrotiger CXX)
set(CMAKE_CXX_STANDARD 20)
add_executable(codegen
  	src/codegen.cpp
   	src/Util.cpp
  	src/Indent.hpp
 	src/Util.hpp
)


