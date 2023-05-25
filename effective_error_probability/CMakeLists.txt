cmake_minimum_required(VERSION 3.24)
# Outputs the CMAKE_BUILD_TYPE
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
# If build type is not set, then use `Debug`
if (CMAKE_BUILD_TYPE STREQUAL "")
  message(STATUS "  Diag: Build type was unspecified, set to Debug")
  set(CMAKE_BUILD_TYPE Debug)
endif()

project(generate_reencoding_errors)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY $<1:${CMAKE_SOURCE_DIR}/build/bin>)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  message(STATUS "  OPENMP: found.")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
if(OPENMP_FOUND)
  file(GLOB SOURCESINFOLDER "src/*.cpp")
  file(GLOB HEADERSINFOLDER "include/*.hpp")
  include_directories(app ${CMAKE_CURRENT_SOURCE_DIR} "./include")
  set(PROJECT_SOURCES
    ${SOURCESINFOLDER}
    ${HEADERSINFOLDER}
  )
  add_executable(main
    ${PROJECT_SOURCES}
  )
else()
    message(FATAL_ERROR "OPENMP not found!")
endif()
