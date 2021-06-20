cmake_minimum_required(VERSION 3.13)

project(newtone-methods)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread")

include_directories("${CMAKE_SOURCE_DIR}/headers")

file(GLOB_RECURSE src "${CMAKE_SOURCE_DIR}/src/*.cpp")

add_executable(newtone-methods "${src}")
