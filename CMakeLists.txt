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


add_executable("asymptotics-1" asymptotics_1.cpp QuantumGraph.cpp)

target_include_directories(asymptotics-1 PUBLIC "${PROJECT_BINARY_DIR}")

target_include_directories(asymptotics-1 PUBLIC "C:/Libs/Mathter/include")

target_include_directories(asymptotics-1 PUBLIC "C:/Libs/boost_1_81_0")

target_link_libraries(asymptotics-1 PUBLIC matplot)


add_executable("asymptotics-2" asymptotics_2.cpp QuantumGraph.cpp)

target_include_directories(asymptotics-2 PUBLIC "${PROJECT_BINARY_DIR}")

target_include_directories(asymptotics-2 PUBLIC "C:/Libs/Mathter/include")

target_include_directories(asymptotics-2 PUBLIC "C:/Libs/boost_1_81_0")


add_executable("determinant" determinant.cpp QuantumGraph.cpp)

target_include_directories(determinant PUBLIC "${PROJECT_BINARY_DIR}")

target_include_directories(determinant PUBLIC "C:/Libs/Mathter/include")

target_include_directories(determinant PUBLIC "C:/Libs/boost_1_81_0")

target_link_libraries(determinant PUBLIC matplot)


add_executable("asymptotics-3" asymptotics_3.cpp QuantumGraph.cpp)

target_include_directories(asymptotics-3 PUBLIC "${PROJECT_BINARY_DIR}")

target_include_directories(asymptotics-3 PUBLIC "C:/Libs/Mathter/include")

target_include_directories(asymptotics-3 PUBLIC "C:/Libs/boost_1_81_0")

target_compile_definitions(asymptotics-3 PUBLIC BOOST_JSON_NO_LIB)


add_executable("eigenvalues" eigenvalues.cpp QuantumGraph.cpp)

target_include_directories(eigenvalues PUBLIC "${PROJECT_BINARY_DIR}")

target_include_directories(eigenvalues PUBLIC "C:/Libs/Mathter/include")

target_include_directories(eigenvalues PUBLIC "C:/Libs/boost_1_81_0")

target_compile_definitions(eigenvalues PUBLIC BOOST_JSON_NO_LIB)


add_executable("dirichlet" dirichlet.cpp QuantumGraph.cpp)

target_include_directories(dirichlet PUBLIC "${PROJECT_BINARY_DIR}")

target_include_directories(dirichlet PUBLIC "C:/Libs/Mathter/include")

target_include_directories(dirichlet PUBLIC "C:/Libs/boost_1_81_0")

target_compile_definitions(dirichlet PUBLIC BOOST_JSON_NO_LIB)