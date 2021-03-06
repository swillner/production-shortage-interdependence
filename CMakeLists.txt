cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(cmake/helpers.cmake)

set_default_build_type(Release)

project(psi)

add_executable(psi src/psi.cpp)

target_include_directories(psi PRIVATE include lib/cpp-library)
target_compile_options(psi PRIVATE "-std=c++11")

set_advanced_cpp_warnings(psi)
set_build_type_specifics(psi)
add_git_version(psi
  FALLBACK_VERSION 1.0.0)

option(PSI_PARALLELIZATION "" ON)
if(PSI_PARALLELIZATION)
  find_package(OpenMP REQUIRED)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include(lib/libmrio/libmrio.cmake)
target_link_libraries(psi libmrio)

add_cpp_tools(psi)
