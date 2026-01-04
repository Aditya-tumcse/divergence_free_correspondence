# Collect all source files for the library
set(SRC_FILES
    io.cpp
    features.cpp
    matching.cpp
    deformation_field.cpp
    utilities.cpp
    numerics.cpp
    run.cpp
)

# Define a static library for the core functionality
add_library(divergence_free_corr STATIC ${SRC_FILES})

# Include directories for the library
target_include_directories(divergence_free_corr PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/external/include
    ${PCL_INCLUDE_DIRS}
    ${Ceres_INCLUDE_DIRS}
)

# Link libraries to the static library
target_link_libraries(divergence_free_corr PUBLIC ${PCL_LIBRARIES} Ceres::ceres)
