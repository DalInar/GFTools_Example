cmake_minimum_required(VERSION 3.5)
project(GFTools_Example)

find_package(ALPSCore COMPONENTS hdf5 accumulators mc params gf REQUIRED)
message(STATUS "ALPSCore includes: " ${ALPSCore_INCLUDE_DIR})
include_directories(${ALPSCore_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(GFTools_Example ${SOURCE_FILES})

target_link_libraries(GFTools_Example
        ${ALPSCore_LIBRARIES}
        )