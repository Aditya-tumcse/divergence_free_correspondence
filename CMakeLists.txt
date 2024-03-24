cmake_minimum_required(VERSION 3.22)
project(divergence_free_corr)

# use C++ 20
set(CMAKE_CXX_STANDARD 20)

# Set include directories for project headers (replace if necessary)
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add the src directory (assuming all source code resides there)
add_subdirectory(src)


