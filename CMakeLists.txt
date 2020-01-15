cmake_minimum_required(VERSION 3.15)
project(ProjectQQMDS_HenonMap)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(SOURCE_FILES main.cpp src/henon.h src/henon.cpp src/stability.cpp src/paramMethod.cpp src/IOmodule.cpp)
add_executable(ProjectQQMDS_HenonMap ${SOURCE_FILES} )

find_package(GSL REQUIRED)
target_link_libraries(ProjectQQMDS_HenonMap GSL::gsl GSL::gslcblas)