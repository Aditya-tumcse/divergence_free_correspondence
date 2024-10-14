# Include directories for PCL
include_directories(${PCL_INCLUDE_DIRS})

# Create an executable
add_executable(divergence_free_corr io.cpp matching.cpp deformation_field.cpp numerics.cpp utilities.cpp run.cpp main.cpp)

# Link PCL libraries to the executable
target_link_libraries(divergence_free_corr PUBLIC ${PCL_LIBRARIES})





