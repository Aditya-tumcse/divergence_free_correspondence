### Divergence Free Correspondence

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

# 1. Computational Efficiency
 - Processing large point clouds within reasonable time constraints using efficient custom implementation of parallelised adjustable farthest point sampling
 - Efficient storage of soft correspondence matrix as a sparse matrix

# 2. Deformation Field Prediction
 - Computation of low frequency velocity basis functions
 - Optimization of velocity field coefficients based on automatic differentiation
 - Tensor representation of velocity basis functions making the computations faster
 - Implementation of custom numerical time integration methods to obtain precise positions of updated point cloud

## Future Improvements
 - [ ] Integration of Open3D library for basic pointcloud processing
 - [ ] Implementation of efficient numerical time integration methods
 - [ ] Improved runtime efficiency of deformation field computation
 - [ ] Enhanced SIMD optimization 
 - [ ] Enhanced GPU utilization

