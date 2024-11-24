#include "deformation_field.hpp"
#include "constants.hpp"
#include "io.hpp"
#include "numerics.hpp"

#include <algorithm>
#include <cmath>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <pthread.h>
#include <thread>

namespace adi {
namespace deformation_field {
Fastor::Tensor<double, 1>
DeformationField::dphidx(const BasisIndices basis_indices,
                         const Eigen::MatrixXd &pts) {

  Eigen::VectorXd dphi_dx;
  dphi_dx.setZero();
  dphi_dx = 0.125 * M_PI * basis_indices.index_1 *
            (M_PI * basis_indices.index_1 * pts.col(0)).array().cos() *
            (M_PI * basis_indices.index_2 * pts.col(1)).array().sin() *
            (M_PI * basis_indices.index_3 * pts.col(2)).array().sin();

  Eigen::Tensor<double, 1> eigen_tensor =
      Eigen::TensorMap<Eigen::Tensor<double, 1>>(dphi_dx.data(),
                                                 dphi_dx.size());
  Fastor::TensorMap<double, 1> fastor_tensor(eigen_tensor.data());

  return fastor_tensor;
}

// template <typename P>
// P DeformationField::dphidx_per_point(const BasisIndices &basis_indices,
//                                      const Eigen::Vector3<P> &point) {
//   P dphidx_per_point;

//   dphidx_per_point = 0.125 * M_PI * basis_indices.index_1 *
//                      ceres::cos(M_PI * basis_indices.index_1 * point[0]) *
//                      ceres::sin(M_PI * basis_indices.index_2 * point[1]) *
//                      ceres::sin(M_PI * basis_indices.index_3 * point[2]);

//   return dphidx_per_point;
// }

Fastor::Tensor<double, 1>
DeformationField::dphidy(const BasisIndices basis_indices,
                         const Eigen::MatrixXd &pts) {
  Eigen::VectorXd dphi_dy;
  dphi_dy.setZero();
  dphi_dy = 0.125 * M_PI * basis_indices.index_2 *
            (M_PI * basis_indices.index_1 * pts.col(0)).array().sin() *
            (M_PI * basis_indices.index_2 * pts.col(1)).array().cos() *
            (M_PI * basis_indices.index_3 * pts.col(2)).array().sin();

  Eigen::Tensor<double, 1> eigen_tensor =
      Eigen::TensorMap<Eigen::Tensor<double, 1>>(dphi_dy.data(),
                                                 dphi_dy.size());
  Fastor::TensorMap<double, 1> fastor_tensor(eigen_tensor.data());

  return fastor_tensor;
}

// template <typename P>
// P DeformationField::dphidy_per_point(const BasisIndices &basis_indices,
//                                      const Eigen::Vector3<P> &point) {
//   P dphidy_per_point;

//   dphidy_per_point = 0.125 * M_PI * basis_indices.index_2 *
//                      ceres::sin(M_PI * basis_indices.index_1 * point[0]) *
//                      ceres::cos(M_PI * basis_indices.index_2 * point[1]) *
//                      ceres::sin(M_PI * basis_indices.index_3 * point[2]);

//   return dphidy_per_point;
// }

Fastor::Tensor<double, 1>
DeformationField::dphidz(const BasisIndices basis_indices,
                         const Eigen::MatrixXd &pts) {
  Eigen::VectorXd dphi_dz;
  dphi_dz.setZero();
  dphi_dz = 0.125 * M_PI * basis_indices.index_3 *
            (M_PI * basis_indices.index_1 * pts.col(0)).array().sin() *
            (M_PI * basis_indices.index_2 * pts.col(1)).array().sin() *
            (M_PI * basis_indices.index_3 * pts.col(2)).array().cos();

  Eigen::Tensor<double, 1> eigen_tensor =
      Eigen::TensorMap<Eigen::Tensor<double, 1>>(dphi_dz.data(),
                                                 dphi_dz.size());
  Fastor::TensorMap<double, 1> fastor_tensor(eigen_tensor.data());

  return fastor_tensor;
}

// template <typename P>
// P DeformationField::dphidz_per_point(const BasisIndices &basis_indices,
//                                      const Eigen::Vector3<P> &point) {
//   P dphidz_per_point;

//   dphidz_per_point = 0.125 * M_PI * basis_indices.index_3 *
//                      ceres::sin(M_PI * basis_indices.index_1 * point[0]) *
//                      ceres::sin(M_PI * basis_indices.index_2 * point[1]) *
//                      ceres::cos(M_PI * basis_indices.index_3 * point[2]);

//   return dphidz_per_point;
// }

Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS, MAX_NUMBER_OF_VELOCITY_BASIS,
               TENSOR_DEPTH>
DeformationField::computeVelocityBasisFunctions(
    const std::vector<BasisIndices> &basis_indices,
    const Eigen::MatrixXd &points) {
  Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS, MAX_NUMBER_OF_VELOCITY_BASIS,
                 TENSOR_DEPTH>
      vel_basis_functions;
  vel_basis_functions.zeros();

  // Compute vectorized partial derivatives of the eigen function
  for (uint32_t i = 0; i < MAX_NUMBER_OF_VELOCITY_BASIS; ++i) {
    vel_basis_functions(Fastor::all, i, 0) = 0.0;
    vel_basis_functions(Fastor::all, i, 1) =
        dphidz(basis_indices.at(i), points);
    vel_basis_functions(Fastor::all, i, 2) =
        -dphidy(basis_indices.at(i), points);
    vel_basis_functions(Fastor::all, i, 3) =
        -vel_basis_functions(Fastor::all, i, 1);
    vel_basis_functions(Fastor::all, i, 4) = 0.0;
    vel_basis_functions(Fastor::all, i, 5) =
        dphidx(basis_indices.at(i), points);
    vel_basis_functions(Fastor::all, i, 6) =
        -vel_basis_functions(Fastor::all, i, 2);
    vel_basis_functions(Fastor::all, i, 7) =
        -vel_basis_functions(Fastor::all, i, 5);
    vel_basis_functions(Fastor::all, i, 8) = 0.0;
  }

  return vel_basis_functions;
}

// template <typename P>
// std::tuple<Eigen::VectorX<P>, Eigen::VectorX<P>, Eigen::VectorX<P>>
// DeformationField::computeVelocityBasisPerPointPerDimension(
//     const Fastor::Tensor<P, 1, MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
//         &vel_basis_per_point) {
//   Fastor::Tensor<P, MAX_NUMBER_OF_VELOCITY_BASIS> vel_basis_per_point_x;
//   Fastor::Tensor<P, MAX_NUMBER_OF_VELOCITY_BASIS> vel_basis_per_point_y;
//   Fastor::Tensor<P, MAX_NUMBER_OF_VELOCITY_BASIS> vel_basis_per_point_z;
//   vel_basis_per_point_x.zeros();
//   vel_basis_per_point_y.zeros();
//   vel_basis_per_point_z.zeros();

//   vel_basis_per_point_x = vel_basis_per_point(Fastor::all, Fastor::all, 0) +
//                           vel_basis_per_point(Fastor::all, Fastor::all, 3) +
//                           vel_basis_per_point(Fastor::all, Fastor::all, 6);
//   vel_basis_per_point_y = vel_basis_per_point(Fastor::all, Fastor::all, 1) +
//                           vel_basis_per_point(Fastor::all, Fastor::all, 4) +
//                           vel_basis_per_point(Fastor::all, Fastor::all, 7);
//   vel_basis_per_point_z = vel_basis_per_point(Fastor::all, Fastor::all, 2) +
//                           vel_basis_per_point(Fastor::all, Fastor::all, 5) +
//                           vel_basis_per_point(Fastor::all, Fastor::all, 8);

//   Eigen::VectorX<P> vel_basis_x = Eigen::Map<Eigen::VectorX<P>>(
//       vel_basis_per_point_x.data(), MAX_NUMBER_OF_VELOCITY_BASIS);
//   Eigen::VectorX<P> vel_basis_y = Eigen::Map<Eigen::VectorX<P>>(
//       vel_basis_per_point_y.data(), MAX_NUMBER_OF_VELOCITY_BASIS);
//   Eigen::VectorX<P> vel_basis_z = Eigen::Map<Eigen::VectorX<P>>(
//       vel_basis_per_point_z.data(), MAX_NUMBER_OF_VELOCITY_BASIS);

//   return std::make_tuple(vel_basis_x, vel_basis_y, vel_basis_z);
// }

// template <typename T>
// std::tuple<Eigen::VectorX<T>, Eigen::VectorX<T>, Eigen::VectorX<T>>
// DeformationField::computeVelocityBasisFunctionsPerPoint(
//     const std::vector<BasisIndices> &basis_indices,
//     const Eigen::Vector3<T> &point) {
//   Fastor::Tensor<T, 1, MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
//       vel_basis_functions_per_point;
//   vel_basis_functions_per_point.zeros();

//   for (uint32_t i = 0; i < MAX_NUMBER_OF_VELOCITY_BASIS; ++i) {
//     vel_basis_functions_per_point(0, i, 0) = 0.0;
//     vel_basis_functions_per_point(0, i, 1) =
//         dphidz_per_point<T>(basis_indices.at(i), point);
//     vel_basis_functions_per_point(0, i, 2) =
//         -dphidy_per_point<T>(basis_indices.at(i), point);
//     vel_basis_functions_per_point(0, i, 3) =
//         -vel_basis_functions_per_point(0, i, 1);
//     vel_basis_functions_per_point(0, i, 4) = 0.0;
//     vel_basis_functions_per_point(0, i, 5) =
//         dphidx_per_point<T>(basis_indices.at(i), point);
//     vel_basis_functions_per_point(0, i, 6) =
//         -vel_basis_functions_per_point(0, i, 2);
//     vel_basis_functions_per_point(0, i, 7) =
//         -vel_basis_functions_per_point(0, i, 5);
//     vel_basis_functions_per_point(0, i, 8) = 0.0;
//   }

//   std::tuple<Eigen::VectorX<T>, Eigen::VectorX<T>, Eigen::VectorX<T>>
//       vel_basis_per_point_per_dimension =
//           computeVelocityBasisPerPointPerDimension<T>(
//               vel_basis_functions_per_point);

//   return vel_basis_per_point_per_dimension;
// }

Eigen::MatrixXd DeformationField::computeVelocityField(
    const Eigen::VectorXd &coeffs_ak,
    const Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
                         MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
        &vel_basis_functions) {
  Eigen::MatrixXd velocity_field_all_pts =
      Eigen::MatrixXd::Zero(NUMBER_OF_SAMPLE_POINTS, 3);

  // Create a Fastor tensor from the Eigen vector of the coefficients of
  // velocity basis functions
  Fastor::Tensor<double, MAX_NUMBER_OF_VELOCITY_BASIS> coeffs_ak_tensor(
      coeffs_ak.data());

  for (uint32_t i = 0; i < TENSOR_DEPTH; i += 3) {
    Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
                   MAX_NUMBER_OF_VELOCITY_BASIS>
        vel_basis_x_component =
            vel_basis_functions(Fastor::all, Fastor::all, i);
    Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
                   MAX_NUMBER_OF_VELOCITY_BASIS>
        intermediate_vel_field_x =
            Fastor::matmul(vel_basis_x_component, coeffs_ak_tensor);

    Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
                   MAX_NUMBER_OF_VELOCITY_BASIS>
        vel_basis_y_component =
            vel_basis_functions(Fastor::all, Fastor::all, i + 1);
    Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
                   MAX_NUMBER_OF_VELOCITY_BASIS>
        intermediate_vel_field_y =
            Fastor::matmul(vel_basis_y_component, coeffs_ak_tensor);

    Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
                   MAX_NUMBER_OF_VELOCITY_BASIS>
        vel_basis_z_component =
            vel_basis_functions(Fastor::all, Fastor::all, i + 2);
    Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
                   MAX_NUMBER_OF_VELOCITY_BASIS>
        intermediate_vel_field_z =
            Fastor::matmul(vel_basis_z_component, coeffs_ak_tensor);

    Eigen::Tensor<double, 2> x_component =
        Eigen::TensorMap<Eigen::Tensor<double, 2>>(
            intermediate_vel_field_x.data(), NUMBER_OF_SAMPLE_POINTS,
            MAX_NUMBER_OF_VELOCITY_BASIS);
    Eigen::Tensor<double, 1> columns_sum_x =
        x_component.sum(Eigen::array<Eigen::Index, 1>({0}));
    Eigen::Map<Eigen::VectorXd> vector_column_sums_x(columns_sum_x.data(),
                                                     columns_sum_x.size());

    Eigen::Tensor<double, 2> y_component =
        Eigen::TensorMap<Eigen::Tensor<double, 2>>(
            intermediate_vel_field_y.data(), NUMBER_OF_SAMPLE_POINTS,
            MAX_NUMBER_OF_VELOCITY_BASIS);
    Eigen::Tensor<double, 1> columns_sum_y =
        y_component.sum(Eigen::array<Eigen::Index, 1>({0}));
    Eigen::Map<Eigen::VectorXd> vector_column_sums_y(columns_sum_y.data(),
                                                     columns_sum_y.size());

    Eigen::Tensor<double, 2> z_component =
        Eigen::TensorMap<Eigen::Tensor<double, 2>>(
            intermediate_vel_field_z.data(), NUMBER_OF_SAMPLE_POINTS,
            MAX_NUMBER_OF_VELOCITY_BASIS);
    Eigen::Tensor<double, 1> columns_sum_z =
        z_component.sum(Eigen::array<Eigen::Index, 1>({0}));
    Eigen::Map<Eigen::VectorXd> vector_column_sums_z(columns_sum_z.data(),
                                                     columns_sum_z.size());

    velocity_field_all_pts.col(0) += vector_column_sums_x;
    velocity_field_all_pts.col(1) += vector_column_sums_y;
    velocity_field_all_pts.col(2) += vector_column_sums_z;
  }

  return velocity_field_all_pts;
}

// template <typename T>
// Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3>
// DeformationField::computeVelocityFieldTemplated(
//     const Eigen::Matrix<T, MAX_NUMBER_OF_VELOCITY_BASIS, 1> &coeffs_ak,
//     const Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
//                          MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
//         &vel_basis_functions) {
//   Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3> velocity_field_all_pts =
//       Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3>::Zero();

//   // Convert coeffs_ak to Fastor tensor compatible with T.
//   Fastor::Tensor<T, MAX_NUMBER_OF_VELOCITY_BASIS> coeffs_ak_tensor(
//       coeffs_ak.data());

//   for (uint32_t i = 0; i < TENSOR_DEPTH; i += 3) {
//     auto vel_basis_x_component =
//         vel_basis_functions(Fastor::all, Fastor::all, i);
//     auto intermediate_vel_field_x =
//         Fastor::matmul(vel_basis_x_component, coeffs_ak_tensor);

//     auto vel_basis_y_component =
//         vel_basis_functions(Fastor::all, Fastor::all, i + 1);
//     auto intermediate_vel_field_y =
//         Fastor::matmul(vel_basis_y_component, coeffs_ak_tensor);

//     auto vel_basis_z_component =
//         vel_basis_functions(Fastor::all, Fastor::all, i + 2);
//     auto intermediate_vel_field_z =
//         Fastor::matmul(vel_basis_z_component, coeffs_ak_tensor);

//     // Map Fastor tensor to Eigen, allowing autodiff types
//     Eigen::TensorMap<Eigen::Tensor<T, 2>> x_component(
//         intermediate_vel_field_x.data(), NUMBER_OF_SAMPLE_POINTS,
//         MAX_NUMBER_OF_VELOCITY_BASIS);
//     Eigen::TensorMap<Eigen::Tensor<T, 2>> y_component(
//         intermediate_vel_field_y.data(), NUMBER_OF_SAMPLE_POINTS,
//         MAX_NUMBER_OF_VELOCITY_BASIS);
//     Eigen::TensorMap<Eigen::Tensor<T, 2>> z_component(
//         intermediate_vel_field_z.data(), NUMBER_OF_SAMPLE_POINTS,
//         MAX_NUMBER_OF_VELOCITY_BASIS);

//     // Sum columns to produce velocity vector components for each point
//     Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 1> vector_column_sums_x =
//         x_component.sum(Eigen::array<Eigen::Index, 1>({0}));
//     Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 1> vector_column_sums_y =
//         y_component.sum(Eigen::array<Eigen::Index, 1>({0}));
//     Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 1> vector_column_sums_z =
//         z_component.sum(Eigen::array<Eigen::Index, 1>({0}));

//     velocity_field_all_pts.col(0) += vector_column_sums_x;
//     velocity_field_all_pts.col(1) += vector_column_sums_y;
//     velocity_field_all_pts.col(2) += vector_column_sums_z;
//   }

//   return velocity_field_all_pts;
// }

} // namespace deformation_field
} // namespace adi