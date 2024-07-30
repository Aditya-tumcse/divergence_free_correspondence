# Inherit project settings from top-level CMakeLists.txt
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(/usr/include/pcl-1.12)
include_directories(/usr/include/eigen3)

# Ensure PCL and Ceres are found and linked
find_package(PCL REQUIRED) 
find_package(Ceres REQUIRED)

# Include PCL directories
include_directories(${PCL_INCLUDE_DIRS})
include_directories(${CERES_INCLUDE_DIRS})

# Create separate executables with their own dependencies
add_executable(src io.cpp matching.cpp deformation_field.cpp numerics.cpp main.cpp)

# Link directories
target_link_libraries(src PUBLIC ${PCL_LIBRARIES} ${CERES_LIBRARIES})




