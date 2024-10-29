# Set source files for the executable
set(SOURCES
    features.cpp
    io.cpp
    matching.cpp
    deformation_field.cpp
    numerics.cpp
    utilities.cpp
    run.cpp
    main.cpp
)

# Create an executable
add_executable(divergence_free_corr ${SOURCES})

# Include directories for PCL and Ceres
target_include_directories(divergence_free_corr PRIVATE ${PCL_INCLUDE_DIRS} ${Ceres_INCLUDE_DIRS} /opt/Fastor)

# Link PCL and Ceres libraries to the executable
target_link_libraries(divergence_free_corr PUBLIC ${PCL_LIBRARIES} ${CERES_LIBRARIES})
