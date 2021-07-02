project(Flimanalysis)
cmake_minimum_required(VERSION 3.19)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()
find_package(GSL REQUIRED)

add_subdirectory(src)