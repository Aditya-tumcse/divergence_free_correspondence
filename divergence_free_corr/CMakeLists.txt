cmake_minimum_required(VERSION 3.16)
project(divergence_free_corr)

# Use C++ 20
set(CMAKE_CXX_STANDARD 20)

# Set include directories for project headers
include_directories(${PROJECT_SOURCE_DIR}/include)

# Find PCL package
find_package(PCL REQUIRED)

# Add the src directory
add_subdirectory(src)
