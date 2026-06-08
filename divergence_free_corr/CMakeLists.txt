# ---------------------------------------------------------------------------
# Core static library
# ---------------------------------------------------------------------------
add_library(divergence_free_corr STATIC
    src/io.cpp
    src/features.cpp
    src/matching.cpp
    src/deformation_field.cpp
    src/utilities.cpp
    src/numerics.cpp
    src/parallel_optimizer.cpp
    src/run.cpp
)

target_include_directories(divergence_free_corr
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/include>
    PRIVATE
        ${PCL_INCLUDE_DIRS}
        ${Ceres_INCLUDE_DIRS}
)

target_link_libraries(divergence_free_corr
    PUBLIC
        ${PCL_LIBRARIES}
        Ceres::ceres
    PRIVATE
        project_options
)

# ---------------------------------------------------------------------------
# Executable
# ---------------------------------------------------------------------------
add_executable(main_executable src/main.cpp)

target_link_libraries(main_executable
    PRIVATE
        divergence_free_corr
        project_options
)
