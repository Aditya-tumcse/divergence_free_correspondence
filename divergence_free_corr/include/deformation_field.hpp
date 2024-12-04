#ifndef DEFORMATION_FIELD_HPP
#define DEFORMATION_FIELD_HPP

#include "constants.hpp"
#include "io.hpp"
#include <Fastor/Fastor.h>

#include <cassert>
#include <ceres/ceres.h>
#include <ceres/jet.h>
#include <eigen3/Eigen/Dense>

namespace adi {
namespace deformation_field {

/**
 * @brief Structure for storing the info regarding the basis indices
 */
struct BasisIndices {
  uint32_t index_1;
  uint32_t index_2;
  uint32_t index_3;
  double eigen_value;
};

template <typename T> double getValue(const T &x) {
  return x; // For non-Jet types (like double)
}

template <typename T, int N> double getValue(const ceres::Jet<T, N> &x) {
  return x.a; // For Jet types
}

class DeformationField {
public:
  DeformationField() {}

  /**
   * @brief Computes the deformation field/velocity vector field for the entire
   * pointcloud
   *
   * @param coeffs_ak
   * @param vel_basis_functions
   *
   * @return An eigen matrix of the velocity field where each column represents
   * each component of the velocity field for each point
   */
  Eigen::MatrixXd
  computeVelocityField(const Eigen::VectorXd &coeffs_ak,
                       const Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
                                            MAX_NUMBER_OF_VELOCITY_BASIS,
                                            TENSOR_DEPTH> &vel_basis_functions);

  // template <typename T>
  // static Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3>
  // computeVelocityFieldTemplated(
  //     const Eigen::Matrix<T, MAX_NUMBER_OF_VELOCITY_BASIS, 1> &coeffs_ak,
  //     const Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
  //                          MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
  //         &vel_basis_functions);

  /**
   * @brief Computes velocity basis functions for the pointcloud
   *
   * @param basis_indices
   * @param points
   *
   * @return Tensor of velocity basis functions for the entire pointcloud and
   * all basis functions
   */
  Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS, MAX_NUMBER_OF_VELOCITY_BASIS,
                 TENSOR_DEPTH>
  computeVelocityBasisFunctions(const std::vector<BasisIndices> &basis_indices,
                                const Eigen::MatrixXd &points);

  template <typename T>
  Eigen::Matrix<T, MAX_NUMBER_OF_VELOCITY_BASIS, 3>
  computeVelocityBasisPerPointPerDimension(
      const Fastor::Tensor<T, 1, MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
          &vel_basis_per_point) {
    Eigen::Matrix<T, MAX_NUMBER_OF_VELOCITY_BASIS, 3> vel_basis_functions;

    for (size_t j = 0; j < MAX_NUMBER_OF_VELOCITY_BASIS; ++j) {
      vel_basis_functions(j, 0) = vel_basis_per_point(0, j, 0) +
                                  vel_basis_per_point(0, j, 3) +
                                  vel_basis_per_point(0, j, 6);

      vel_basis_functions(j, 1) = vel_basis_per_point(0, j, 1) +
                                  vel_basis_per_point(0, j, 4) +
                                  vel_basis_per_point(0, j, 7);

      vel_basis_functions(j, 2) = vel_basis_per_point(0, j, 2) +
                                  vel_basis_per_point(0, j, 5) +
                                  vel_basis_per_point(0, j, 8);
    }
    return vel_basis_functions;
  }

  template <typename T>
  Eigen::Matrix<T, MAX_NUMBER_OF_VELOCITY_BASIS, 3>
  computeVelocityBasisFunctionsPerPoint(
      const std::vector<BasisIndices> &basis_indices,
      const Eigen::Vector3<T> &point) {
    Fastor::Tensor<T, 1, MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
        vel_basis_functions_per_point;

    for (uint32_t i = 0; i < MAX_NUMBER_OF_VELOCITY_BASIS; ++i) {
      vel_basis_functions_per_point(0, i, 0) = T(0.0);
      vel_basis_functions_per_point(0, i, 1) =
          dphidz_per_point<T>(basis_indices.at(i), point);
      vel_basis_functions_per_point(0, i, 2) =
          -dphidy_per_point<T>(basis_indices.at(i), point);
      vel_basis_functions_per_point(0, i, 3) =
          -vel_basis_functions_per_point(0, i, 1);
      vel_basis_functions_per_point(0, i, 4) = T(0.0);
      vel_basis_functions_per_point(0, i, 5) =
          dphidx_per_point<T>(basis_indices.at(i), point);
      vel_basis_functions_per_point(0, i, 6) =
          -vel_basis_functions_per_point(0, i, 2);
      vel_basis_functions_per_point(0, i, 7) =
          -vel_basis_functions_per_point(0, i, 5);
      vel_basis_functions_per_point(0, i, 8) = T(0.0);
    }

    Eigen::Matrix<T, MAX_NUMBER_OF_VELOCITY_BASIS, 3>
        vel_basis_per_point_per_dimension =
            computeVelocityBasisPerPointPerDimension(
                vel_basis_functions_per_point);

    return vel_basis_per_point_per_dimension;
  }

  ~DeformationField() = default;

private:
  /**
   * @brief Function to compute partial derivative of scalar potential field
   * with respect to x at every point in the domain
   *
   * @param basis_indices
   * @param pts
   *
   * @return Tensor of 1 dimesion that repreents partial derivative of the eigen
   * function with respect to x
   */
  Fastor::Tensor<double, 1> dphidx(const BasisIndices basis_indices,
                                   const Eigen::MatrixXd &pts);
  template <typename T>
  T dphidx_per_point(const BasisIndices &basis_indices,
                     const Eigen::Vector3<T> &point) {
    T dphidx_per_point;

    dphidx_per_point =
        T(0.125) * T(M_PI) * T(basis_indices.index_1) *
        ceres::cos(T(M_PI) * T(basis_indices.index_1) * T(point[0])) *
        ceres::sin(T(M_PI) * T(basis_indices.index_2) * T(point[1])) *
        ceres::sin(T(M_PI) * T(basis_indices.index_3) * T(point[2]));

    return dphidx_per_point;
  }

  /**
   * @brief Function to compute partial derivative of scalar potential field
   * with respect to x at every point in the domain
   *
   * @param basis_indices
   * @param pts
   *
   * @return Tensor of 1 dimesion that repreents partial derivative of the eigen
   * function with respect to y
   */
  Fastor::Tensor<double, 1> dphidy(const BasisIndices basis_indices,
                                   const Eigen::MatrixXd &pts);

  template <typename T>
  T dphidy_per_point(const BasisIndices &basis_indices,
                     const Eigen::Vector3<T> &point) {
    T dphidy_per_point;

    dphidy_per_point =
        T(0.125) * T(M_PI) * T(basis_indices.index_2) *
        ceres::sin(T(M_PI) * T(basis_indices.index_1) * T(point[0])) *
        ceres::cos(T(M_PI) * T(basis_indices.index_2) * T(point[1])) *
        ceres::sin(T(M_PI) * T(basis_indices.index_3) * T(point[2]));

    return dphidy_per_point;
  }

  /**
   * @brief Function to compute partial derivative of scalar potential field
   * with respect to x at every point in the domain
   *
   * @param basis_indices
   * @param pts
   *
   * @return Tensor of 1 dimesion that repreents partial derivative of the eigen
   * function with respect to z
   */
  Fastor::Tensor<double, 1> dphidz(const BasisIndices basis_indices,
                                   const Eigen::MatrixXd &pts);

  template <typename T>
  T dphidz_per_point(const BasisIndices &basis_indices,
                     const Eigen::Vector3<T> &point) {
    T dphidz_per_point;

    dphidz_per_point =
        T(0.125) * T(M_PI) * T(basis_indices.index_3) *
        ceres::sin(T(M_PI) * T(basis_indices.index_1) * T(point[0])) *
        ceres::sin(T(M_PI) * T(basis_indices.index_2) * T(point[1])) *
        ceres::cos(T(M_PI) * T(basis_indices.index_3) * T(point[2]));

    return dphidz_per_point;
  }
};
} // namespace deformation_field
} // namespace adi
#endif
