cmake_minimum_required(VERSION 3.24)
# Outputs the CMAKE_BUILD_TYPE
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
# If build type is not set, then use `Debug`
if (CMAKE_BUILD_TYPE STREQUAL "")
  message(STATUS "  Diag: Build type was unspecified, set to Debug")
  set(CMAKE_BUILD_TYPE Debug)
endif()

project(secret_key_rate)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(Boost_DETAILED_FAILURE_MSG ON)
if (WIN32)
    add_definitions(-DBOOST_ALL_NO_LIB)
    set(BOOST_ROOT "C:/local/boost_1_82_0")
    set(BOOST_INCLUDE_DIR "C:/local/boost_1_82_0")
    set(BOOST_LIBRARY_DIR "C:/local/boost_1_82_0/lib64-msvc-14.3")
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON) # https://stackoverflow.com/questions/1941443/cmake-linking-against-shared-library-on-windows-error-about-not-finding-lib-fi
endif ()

set(Boost_USE_MULTITHREADED TRUE)
set(USE_STATIC_BOOST TRUE)
set(Boost_USE_STATIC_LIBS TRUE)
find_package(Boost 1.82.0 COMPONENTS program_options REQUIRED)

if (APPLE)
  find_package(fmt) # brew install fmt
  if(FMT_FOUND)
    message(FATAL_ERROR "fmt not found!")
  endif()
endif()

if(Boost_FOUND)
  file(GLOB SOURCESINFOLDER "src/*.cpp")
  file(GLOB HEADERSINFOLDER "include/*.hpp")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${BOOST_INCLUDE_DIR} "./include" "../effective_error_probability/include")
  set(PROJECT_SOURCES
    ../effective_error_probability/src/AuxMath.cpp
    ${SOURCESINFOLDER}
    ${HEADERSINFOLDER}
  )
  add_executable(main
    ${PROJECT_SOURCES}
  )
  if (WIN32)
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
      target_link_libraries(main optimized ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE})
    else()
      target_link_libraries(main debug ${Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG})
    endif()
  else ()
      if (CMAKE_BUILD_TYPE STREQUAL "Release")
      target_link_libraries(main optimized ${Boost_LIBRARIES})
    else()
      target_link_libraries(main debug ${Boost_LIBRARIES})
    endif()
  endif()
else()
  message(FATAL_ERROR "Boost not found!")
endif()

