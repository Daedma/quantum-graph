cmake_minimum_required(VERSION 3.26.4)

project(quantum-graph VERSION 1.0.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT DEFINED QGRAPH_DEPS_DIRECTORY)
	set(QGRAPH_DEPS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/deps)
endif()

find_package(Boost REQUIRED COMPONENTS odeint)

add_subdirectory(${QGRAPH_DEPS_DIRECTORY}/Mathter)

add_library(qgraph QuantumGraph.cpp)

target_include_directories(qgraph PUBLIC 
	${CMAKE_CURRENT_SOURCE_DIR} 
	${Boost_INCLUDE_DIRS} 
	${QGRAPH_DEPS_DIRECTORY}/Mathter/include)