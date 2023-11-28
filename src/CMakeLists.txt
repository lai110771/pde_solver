cmake_minimum_required(VERSION 3.1)
project(PDE-Solver VERSION 0.1.0 LANGUAGES CXX)
find_package (Eigen3 3.3 REQUIRED NO_MUDULE)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/config.txt
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_executable(main main.cpp Coordinate.cpp Solver.cpp)
target_link_libraries (main Eigen3::Eigen)
set (CMAKE_CXX_FLAGS "-O3")
