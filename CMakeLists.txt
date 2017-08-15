cmake_minimum_required(VERSION 2.8)
if(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.8)
  cmake_policy(SET CMP0069 NEW) # for INTERPROCEDURAL_OPTIMIZATION
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_SOURCE_DIR}/cmake/find-modules)

project(acclimate)

set(PSI_VERSION 1.0.0)
if(EXISTS "${CMAKE_SOURCE_DIR}/.git" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/.git")
  find_program(HAVE_GIT git)
  if(HAVE_GIT)
    execute_process(COMMAND git describe --tags --dirty --always OUTPUT_VARIABLE GIT_OUTPUT OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX REPLACE "^v([0-9]+\\.[0-9]+)\\.(0-)?([0-9]*)((-.+)?)$" "\\1.\\3\\4" PSI_VERSION "${GIT_OUTPUT}")
  endif()
  mark_as_advanced(HAVE_GIT)
endif()
message(STATUS "Program version: ${PSI_VERSION}")

add_executable(psi src/psi.cpp)

target_include_directories(psi PRIVATE include lib/cpp-library)
target_compile_options(psi PRIVATE "-std=c++11")
target_compile_definitions(psi PRIVATE PSI_VERSION="${PSI_VERSION}")

if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" OR CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
  if(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.8)
    message(STATUS "Enabling interprocedural optimization")
    set_property(TARGET psi PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
  target_compile_definitions(psi PUBLIC NDEBUG)
else()
  target_compile_definitions(psi PRIVATE DEBUG)
endif()

#include NetCDF C library
find_package(NETCDF REQUIRED)
message(STATUS "NetCDF include directory: ${NETCDF_INCLUDE_DIR}")
message(STATUS "NetCDF library: ${NETCDF_LIBRARY}")
target_link_libraries(psi netcdf)

#include NetCDF C++ library
find_package(NETCDF_CPP4 REQUIRED)
message(STATUS "NetCDF_c++4 include directory: ${NETCDF_CPP4_INCLUDE_DIR}")
message(STATUS "NetCDF_c++4 library: ${NETCDF_CPP4_LIBRARY}")
target_link_libraries(psi netcdf_c++4)

include(lib/libmrio/libmrio.cmake)
target_link_libraries(psi libmrio)
