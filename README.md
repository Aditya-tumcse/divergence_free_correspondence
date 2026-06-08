# Divergence Free Correspondence

## Introduction

DVF provides a robust solution for manipulating and deforming 3D point cloud data with precision and efficiency. The library implements advanced velocity field-based deformation model optimized through a parallel correspondence computation pipeline.

## Main Features

- Computation of SHOT features to obtain initial correspondences
- Fast parallel optimization of deformation fields
- Efficient sparse matrix correspondence computation
- Deformation field prediction using the concept of divergence free velocity fields
- Compatible with Ceres Solver for non-linear optimization
- Memory efficient processing of large point clouds
- Multi-threaded implementation for maximum performance

## Technical Challenges and Solutions
This section documents the major technical challenges encountered and the approaches used to solve them

### 1. Computational Efficiency
 - Processing large point clouds within reasonable time constraints using efficient custom implementation of parallelised adjustable farthest point sampling
 - Efficient storage of soft correspondence matrix as a sparse matrix

### 2. Deformation Field Prediction
 - Computation of low frequency velocity basis functions
 - Optimization of velocity field coefficients based on automatic differentiation
 - Tensor representation of velocity basis functions making the computations faster
 - Implementation of custom numerical time integration methods to obtain precise positions of updated point cloud

## Project Structure

```
divergence_free_correspondence/
├ divergence_free_corr/
| ├ external/                   # External dependencies folder
| | ├ Fastor/                   # Fastor repo header files
| | └ clue.hpp                  # Logging repo header files
| ├ include/                    # Folder containing the header files
| | ├ constants.hpp             # Contains the constant expressions used in the project
| | ├ deformation_field.hpp     # Deformation field related computations 
| | ├ features.hpp              # SHOT features extraction
| | ├ io.hpp                    # Reading and processing of point clouds
| | ├ matching.hpp              # Descriptor matching computation
| | ├ numerics.hpp              # Numerical time integration schemes
| | ├ parallel_optimizer.hpp    # Parallel implementation of optimizer
| | ├ run.hpp                   # Header for the run function
| | └ utilities.hpp             # Utility components
| ├ src/                        # Folder containing the corresponding cpp files
| | ├ CMakeLists.txt
| | ├ deformation_field.cpp
| | ├ features.cpp
| | ├ io.cpp
| | ├ main.cpp
| | ├ matching.cpp
| | ├ numerics.cpp
| | ├ parallel_optimizer.cpp
| | ├ run.cpp
| | └ utilities.cpp
| └ CMakeLists.txt              # Main CMake configuration file
├ data/                         # Folder containing the input point clouds
| ├ source_cloud.ply            # Source point cloud
| └ target_cloud.ply            # Target point cloud
├ docker/                       # Docker folder
| ├ build_dependencies.sh       # Bash file to build docker dependencies for the repo 
| ├ Dockerfile_dependencies     # Docker file containing the dependencies to be built to run the repo
| └ Dockerfile_develop          # Docker file needed to build the dev container
├ .devcontainer                 # Dev container configuration folder
├ .gitignore                    # Git ignore file
├ .clang-format                 # Clang format file             
└ README.md                     # Project overview and documentation

```

## Requirements
 - CMake >= 3.16
 - Python >= 3.8
 - [Conan](https://conan.io/) >= 2.0 (`pip install conan`)

## Installation

### With Dev Container (recommended)

Open the repo in VS Code and reopen in the dev container. The `postCreateCommand` will automatically:
1. Detect your compiler profile
2. Fetch and build all C++ dependencies via Conan
3. Generate the CMake toolchain file

### Local build

```bash
# Install Conan if not already installed
pip install conan

# Detect compiler profile (once per machine)
conan profile detect --force

# Fetch/build dependencies (output goes into build/)
conan install . --output-folder=build --build=missing \
    -s build_type=Release -s compiler.cppstd=20

# Configure and build
cmake -S divergence_free_corr -B build \
    -DCMAKE_TOOLCHAIN_FILE=build/conan_toolchain.cmake \
    -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

> **Note:** The first `conan install` will build PCL and Ceres from source if pre-built
> binaries are not available for your compiler profile. This can take 20–40 minutes.
> Subsequent runs use the Conan cache and are fast.

To check available PCL versions on ConanCenter:
```bash
conan search pcl -r conancenter
```

## Usage
Ensure there is a data folder as shown in the folder structure above.
Run the executable using the following command

```bash
./build/main_executable
```

## Future Improvements
 - [ ] Integration of Open3D library for basic pointcloud processing
 - [ ] Implementation of efficient numerical time integration methods
 - [ ] Improved runtime efficiency of deformation field computation
 - [ ] Build python bindings for the repo
 - [ ] Enhanced SIMD optimization 
 - [ ] Enhanced GPU utilization

## About

This project was created by Aditya Sai Srinivas.

