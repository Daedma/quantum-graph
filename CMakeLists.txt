cmake_minimum_required(VERSION 3.26.4)

project("Quantum graphs" VERSION 0.0.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(matplotplusplus)

add_executable("qgraphs" main.cpp)

target_include_directories(qgraphs PUBLIC "${PROJECT_BINARY_DIR}")

target_include_directories(qgraphs PUBLIC "C:/Libs/boost_1_81_0")

target_link_libraries(qgraphs PUBLIC matplot)