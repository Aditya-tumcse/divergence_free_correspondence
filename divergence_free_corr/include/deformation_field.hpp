#ifndef DEFORMATION_FIELD_HPP
#define DEFORMATION_FIELD_HPP

#include "constants.hpp"
#include "io.hpp"
#include <Fastor/Fastor.h>

#include <cassert>
#include <ceres/ceres.h>
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
  std::tuple<Eigen::VectorX<T>, Eigen::VectorX<T>, Eigen::VectorX<T>>
  computeVelocityBasisPerPointPerDimension(
      const Fastor::Tensor<T, 1, MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
          &vel_basis_per_point) {
    // Fastor::Tensor<T, MAX_NUMBER_OF_VELOCITY_BASIS> vel_basis_per_point_x;
    // Fastor::Tensor<T, MAX_NUMBER_OF_VELOCITY_BASIS> vel_basis_per_point_y;
    // Fastor::Tensor<T, MAX_NUMBER_OF_VELOCITY_BASIS> vel_basis_per_point_z;
    // vel_basis_per_point_x.zeros();
    // vel_basis_per_point_y.zeros();
    // vel_basis_per_point_z.zeros();
    Eigen::VectorX<T> vel_basis_x;
    vel_basis_x.setZero();
    Eigen::VectorX<T> vel_basis_y;
    vel_basis_y.setZero();
    Eigen::VectorX<T> vel_basis_z;
    vel_basis_z.setZero();

    // vel_basis_per_point_x = vel_basis_per_point(Fastor::all, Fastor::all, 0)
    // +
    //                         vel_basis_per_point(Fastor::all, Fastor::all, 3)
    //                         + vel_basis_per_point(Fastor::all, Fastor::all,
    //                         6);
    // vel_basis_per_point_y = vel_basis_per_point(Fastor::all, Fastor::all, 1)
    // +
    //                         vel_basis_per_point(Fastor::all, Fastor::all, 4)
    //                         + vel_basis_per_point(Fastor::all, Fastor::all,
    //                         7);
    // vel_basis_per_point_z = vel_basis_per_point(Fastor::all, Fastor::all, 2)
    // +
    //                         vel_basis_per_point(Fastor::all, Fastor::all, 5)
    //                         + vel_basis_per_point(Fastor::all, Fastor::all,
    //                         8);

    // Eigen::VectorX<T> vel_basis_x = Eigen::Map<Eigen::VectorX<T>>(
    //     vel_basis_per_point_x.data(), MAX_NUMBER_OF_VELOCITY_BASIS);
    // Eigen::VectorX<T> vel_basis_y = Eigen::Map<Eigen::VectorX<T>>(
    //     vel_basis_per_point_y.data(), MAX_NUMBER_OF_VELOCITY_BASIS);
    // Eigen::VectorX<T> vel_basis_z = Eigen::Map<Eigen::VectorX<T>>(
    //     vel_basis_per_point_z.data(), MAX_NUMBER_OF_VELOCITY_BASIS);
    for (size_t j = 0; j < MAX_NUMBER_OF_VELOCITY_BASIS; ++j) {
      vel_basis_x(j) = vel_basis_per_point(0, j, 0) +
                       vel_basis_per_point(0, j, 3) +
                       vel_basis_per_point(0, j, 6);

      vel_basis_y(j) = vel_basis_per_point(0, j, 1) +
                       vel_basis_per_point(0, j, 4) +
                       vel_basis_per_point(0, j, 7);

      vel_basis_z(j) = vel_basis_per_point(0, j, 2) +
                       vel_basis_per_point(0, j, 5) +
                       vel_basis_per_point(0, j, 8);
    }

    return std::make_tuple(vel_basis_x, vel_basis_y, vel_basis_z);
  }

  template <typename T>
  std::tuple<Eigen::VectorX<T>, Eigen::VectorX<T>, Eigen::VectorX<T>>
  computeVelocityBasisFunctionsPerPoint(
      const std::vector<BasisIndices> &basis_indices,
      const Eigen::Vector3<T> &point) {
    Fastor::Tensor<T, 1, MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
        vel_basis_functions_per_point;
    // vel_basis_functions_per_point.zeros();

    for (uint32_t i = 0; i < MAX_NUMBER_OF_VELOCITY_BASIS; ++i) {
      vel_basis_functions_per_point(0, i, 0) = T(0);
      vel_basis_functions_per_point(0, i, 1) =
          dphidz_per_point<T>(basis_indices.at(i), point);
      vel_basis_functions_per_point(0, i, 2) =
          -dphidy_per_point<T>(basis_indices.at(i), point);
      vel_basis_functions_per_point(0, i, 3) =
          vel_basis_functions_per_point(0, i, 4) = T(0);
      -vel_basis_functions_per_point(0, i, 1);
      vel_basis_functions_per_point(0, i, 5) =
          dphidx_per_point<T>(basis_indices.at(i), point);
      vel_basis_functions_per_point(0, i, 6) =
          -vel_basis_functions_per_point(0, i, 2);
      vel_basis_functions_per_point(0, i, 7) =
          -vel_basis_functions_per_point(0, i, 5);
      vel_basis_functions_per_point(0, i, 8) = T(0);
    }

    std::tuple<Eigen::VectorX<T>, Eigen::VectorX<T>, Eigen::VectorX<T>>
        vel_basis_per_point_per_dimension =
            computeVelocityBasisPerPointPerDimension<T>(
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

  template <typename P>
  P dphidx_per_point(const BasisIndices &basis_indices,
                     const Eigen::Vector3<P> &point) {
    P dphidx_per_point;

    dphidx_per_point = 0.125 * M_PI * basis_indices.index_1 *
                       ceres::cos(M_PI * basis_indices.index_1 * point[0]) *
                       ceres::sin(M_PI * basis_indices.index_2 * point[1]) *
                       ceres::sin(M_PI * basis_indices.index_3 * point[2]);

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

  template <typename P>
  P dphidy_per_point(const BasisIndices &basis_indices,
                     const Eigen::Vector3<P> &point) {
    P dphidy_per_point;

    dphidy_per_point = 0.125 * M_PI * basis_indices.index_2 *
                       ceres::sin(M_PI * basis_indices.index_1 * point[0]) *
                       ceres::cos(M_PI * basis_indices.index_2 * point[1]) *
                       ceres::sin(M_PI * basis_indices.index_3 * point[2]);

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

  template <typename P>
  P dphidz_per_point(const BasisIndices &basis_indices,
                     const Eigen::Vector3<P> &point) {
    P dphidz_per_point;

    dphidz_per_point = 0.125 * M_PI * basis_indices.index_3 *
                       ceres::sin(M_PI * basis_indices.index_1 * point[0]) *
                       ceres::sin(M_PI * basis_indices.index_2 * point[1]) *
                       ceres::cos(M_PI * basis_indices.index_3 * point[2]);

    return dphidz_per_point;
  }
};
} // namespace deformation_field
} // namespace adi
#endif
