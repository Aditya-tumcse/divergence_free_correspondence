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
} // namespace deformation_field
} // namespace adi