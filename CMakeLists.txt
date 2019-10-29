cmake_minimum_required(VERSION 3.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(cmake/helpers.cmake)
include(cmake/libraries.cmake)

set_default_build_type(Release)

project(psi)

add_executable(psi src/psi.cpp)
target_include_directories(psi PRIVATE include lib/cpp-library lib/nvector)
set_property(TARGET psi PROPERTY CXX_STANDARD 11)

add_system_paths()
set_advanced_options()
set_advanced_cpp_warnings(psi)
set_build_type_specifics(psi)
add_git_version(psi)

option(PSI_PARALLELIZATION "" ON)
if(PSI_PARALLELIZATION)
  find_package(OpenMP REQUIRED)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include(lib/libmrio/libmrio.cmake)
target_link_libraries(psi PRIVATE libmrio)
mark_as_advanced(LIBMRIO_WITH_NETCDF)

add_cpp_tools(psi)
