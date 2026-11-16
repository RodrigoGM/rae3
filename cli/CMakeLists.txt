cmake_minimum_required(VERSION 3.16)
project(rae LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
# -ffp-contract=off keeps floating-point results identical across compilers and with the R build.
add_compile_options(-ffp-contract=off)
add_executable(rae main.cpp)
# The header-only core lives inside the R package so both front ends compile the same code.
target_include_directories(rae PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../rae3/inst/include)
install(TARGETS rae RUNTIME DESTINATION bin)
