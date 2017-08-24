cmake_minimum_required(VERSION 3.4)

project(rpyfmm CXX) 

set(CMAKE_MODULE_PATH ${CMAKE_MOUDLE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3") 

include(DownloadProject) 

download_project(
  PROJ dashmm
  GIT_REPOSITORY https://github.com/jacksondebuhr/dashmm.git
  GIT_TAG dev-v2
  ${UPDATE_DISCONNECTED_IF_AVAILABLE})
add_subdirectory(${dashmm_SOURCE_DIR} ${dashmm_BINARY_DIR})

find_package(PkgConfig)
pkg_search_module(HPX REQUIRED hpx)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/example)

#add_subdirectory(src)
#add_subdirectory(include)

