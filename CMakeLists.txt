cmake_minimum_required(VERSION 3.26.4)

project("Quantum graphs" VERSION 0.3.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(matplotplusplus-1.2.1)

add_executable("error_test" error_test.cpp QuantumGraph.cpp)

target_include_directories(error_test PUBLIC "${PROJECT_BINARY_DIR}")

target_include_directories(error_test PUBLIC "C:/Libs/Mathter/include")

target_include_directories(error_test PUBLIC "C:/Libs/boost_1_81_0")

add_executable("eigenfunctions" eigenfunctions.cpp QuantumGraph.cpp)

target_include_directories(eigenfunctions PUBLIC "${PROJECT_BINARY_DIR}")

target_include_directories(eigenfunctions PUBLIC "C:/Libs/Mathter/include")

target_include_directories(eigenfunctions PUBLIC "C:/Libs/boost_1_81_0")

target_link_libraries(eigenfunctions PUBLIC matplot)

add_executable("sequences" sequences.cpp QuantumGraph.cpp)

target_include_directories(sequences PUBLIC "${PROJECT_BINARY_DIR}")

target_include_directories(sequences PUBLIC "C:/Libs/Mathter/include")

target_include_directories(sequences PUBLIC "C:/Libs/boost_1_81_0")
