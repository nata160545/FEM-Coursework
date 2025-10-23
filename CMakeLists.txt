cmake_minimum_required(VERSION 3.10)
project(FEM_Solver)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add executable
add_executable(fem_solver main.cpp)

# Include Eigen headers
target_include_directories(fem_solver PRIVATE 
    "${CMAKE_SOURCE_DIR}/include"
)

# Windows specific settings
if(WIN32)
    target_compile_definitions(fem_solver PRIVATE _USE_MATH_DEFINES)
endif()

message(STATUS "Include path: ${CMAKE_SOURCE_DIR}/include")