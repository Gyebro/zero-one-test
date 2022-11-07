cmake_minimum_required(VERSION 3.20)
project(zero_one_test)

set(CMAKE_CXX_STANDARD 14)

set(SOURCES
        src/gnuplot_i.cpp
        src/vec2.cpp
        src/system.cpp
        src/microchaos.cpp
        src/main.cpp)

add_executable(zero_one_test ${SOURCES})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(zero_one_test PUBLIC OpenMP::OpenMP_CXX)
endif()
