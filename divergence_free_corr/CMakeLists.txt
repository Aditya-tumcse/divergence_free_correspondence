cmake_minimum_required(VERSION 3.16)
project(divergence_free_corr)

# Use C++ 20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set include directories for project headers
include_directories(${PROJECT_SOURCE_DIR}/include)

# Find required packages
find_package(PCL REQUIRED)
find_package(Ceres REQUIRED)

# Add the `src` subdirectory
add_subdirectory(src)

# Create the main executable
add_executable(main_executable src/main.cpp)

# Link PCL and other dependencies to the main executable
target_include_directories(main_executable PRIVATE ${PCL_INCLUDE_DIRS} ${Ceres_INCLUDE_DIRS})
target_link_libraries(main_executable PRIVATE divergence_free_corr ${PCL_LIBRARIES} Ceres::ceres)



