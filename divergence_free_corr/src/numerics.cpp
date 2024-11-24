#include "numerics.hpp"
#include "../external/include/clue.hpp"
#include "utilities.hpp"

#include <Fastor/Fastor.h>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <random>

namespace adi {
namespace numerics {
const std::vector<adi::deformation_field::BasisIndices>
GenerateBasisIndices(const uint32_t &max_number_of_velocity_basis) {
  std::vector<adi::deformation_field::BasisIndices> basis_indices;
  double pi_squared = std::pow(M_PI, 2);
  for (uint32_t i = 1; i <= max_number_of_velocity_basis; ++i) {
    double i_squared = std::pow(i, 2);
    for (uint32_t j = 1; j <= max_number_of_velocity_basis; ++j) {
      double j_squared = std::pow(j, 2);
      for (uint32_t k = 1; k <= max_number_of_velocity_basis; ++k) {
        adi::deformation_field::BasisIndices base_index;
        double k_squared = std::pow(k, 2);
        base_index.index_1 = i;
        base_index.index_2 = j;
        base_index.index_3 = k;
        double eigen_val =
            std::pow((pi_squared * (i_squared + j_squared + k_squared)), -1.5);
        base_index.eigen_value = eigen_val;

        basis_indices.push_back(base_index);
      }
    }
  }

  // Sort the basis indices in ascending order based on eigen_value
  std::sort(basis_indices.begin(), basis_indices.end(),
            [](const adi::deformation_field::BasisIndices &a,
               const adi::deformation_field::BasisIndices &b) {
              return a.eigen_value < b.eigen_value;
            });

  return basis_indices;
}

const Eigen::MatrixXd computePrecisionMatrix(
    const std::vector<adi::deformation_field::BasisIndices> &basis_indices) {
  Eigen::MatrixXd precision_matrix =
      Eigen::MatrixXd::Identity(basis_indices.size(), basis_indices.size());

  Eigen::VectorXd inverse_eigen_vals =
      Eigen::VectorXd::Zero(basis_indices.size());

  for (uint32_t i = 0; i < basis_indices.size(); ++i) {
    inverse_eigen_vals[i] = 1.0 / (basis_indices.at(i).eigen_value);
  }

  return (precision_matrix * inverse_eigen_vals);
}

Eigen::MatrixXd RungeKutta2Integration(
    const std::vector<adi::deformation_field::BasisIndices> &basis_indices,
    Eigen::MatrixXd src_point_cloud, const Eigen::VectorXd &coeffs_ak,
    const Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
                         MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
        &vel_basis_functions,
    const uint32_t num_time_steps) {
  Eigen::MatrixXd updated_pts = src_point_cloud;

  const double dt = 1.0 / num_time_steps;

  for (uint32_t t = 0; t < num_time_steps; ++t) {
    adi::deformation_field::DeformationField df;
    Eigen::MatrixXd intermediate_vel_field =
        df.computeVelocityField(coeffs_ak, vel_basis_functions);

    Eigen::MatrixXd midpoint =
        Eigen::MatrixXd::Zero(NUMBER_OF_SAMPLE_POINTS, 3);

    midpoint.col(0) =
        updated_pts.col(0) + 0.5 * dt * intermediate_vel_field.col(0);
    midpoint.col(1) =
        updated_pts.col(1) + 0.5 * dt * intermediate_vel_field.col(1);
    midpoint.col(2) =
        updated_pts.col(2) + 0.5 * dt * intermediate_vel_field.col(2);

    auto updated_vel_basis_functions =
        df.computeVelocityBasisFunctions(basis_indices, midpoint);

    Eigen::MatrixXd updated_vel_field =
        df.computeVelocityField(coeffs_ak, updated_vel_basis_functions);

    updated_pts.col(0) += dt * updated_vel_field.col(0);
    updated_pts.col(1) += dt * updated_vel_field.col(1);
    updated_pts.col(2) += dt * updated_vel_field.col(2);
  }

  return updated_pts;
}

const double EucledianDistance(const Eigen::Vector3d &point_1,
                               const Eigen::Vector3d &point_2) {
  return (point_1 - point_2).norm();
}

template <typename T>
Eigen::Vector3<T> PositionIncrementor<T>::computeVelocityField(
    const Eigen::VectorX<T> &vel_basis_x, const Eigen::VectorX<T> &vel_basis_y,
    const Eigen::VectorX<T> &vel_basis_z, const Eigen::VectorX<T> &coeffs) {
  Eigen::Vector3<T> vel_field;
  vel_field.setZero();

  vel_field[0] = T(vel_basis_x.transpose() * coeffs);
  vel_field[1] = T(vel_basis_y.transpose() * coeffs);
  vel_field[2] = T(vel_basis_z.transpose() * coeffs);

  return vel_field;
}

template <typename T>
void PositionIncrementor<T>::incrementPosition(
    const Eigen::VectorX<T> &vel_basis_x, const Eigen::VectorX<T> &vel_basis_y,
    const Eigen::VectorX<T> &vel_basis_z, const Eigen::VectorX<T> &coeffs,
    const std::vector<adi::deformation_field::BasisIndices> &basis_indices,
    const double num_time_steps) {

  const T dt = T(1.0) / T(num_time_steps);

  for (uint32_t t = 0; t < num_time_steps; ++t) {

    Eigen::Vector3<T> intermediate_vel_field = this->computeVelocityField(
        vel_basis_x, vel_basis_y, vel_basis_z, coeffs);

    (*s_updated_pt)[0] += T(0.5) * dt * intermediate_vel_field[0];
    (*s_updated_pt)[1] += T(0.5) * dt * intermediate_vel_field[1];
    (*s_updated_pt)[2] += T(0.5) * dt * intermediate_vel_field[2];

    auto [vel_basis_function_x, vel_basis_function_y, vel_basis_function_z] =
        s_df.computeVelocityBasisFunctionsPerPoint<T>(basis_indices,
                                                      *s_updated_pt);

    Eigen::Vector3<T> final_vel_field =
        this->computeVelocityField(vel_basis_function_x, vel_basis_function_y,
                                   vel_basis_function_z, coeffs);

    (*s_updated_pt)[0] += dt * final_vel_field[0];
    (*s_updated_pt)[1] += dt * final_vel_field[1];
    (*s_updated_pt)[2] += dt * final_vel_field[2];
  }
}

struct MatchingCostFunctor {
  MatchingCostFunctor(
      const Eigen::Vector3d &source_point, const Eigen::Vector3d &target_point,
      const Eigen::VectorXd &vel_basis_functions_x,
      const Eigen::VectorXd &vel_basis_functions_y,
      const Eigen::VectorXd &vel_basis_functions_z, const double weight,
      const std::vector<adi::deformation_field::BasisIndices> &base_indices)
      : s_target_point(target_point), s_source_point(source_point),
        s_vel_basis_function_x(vel_basis_functions_x),
        s_vel_basis_function_y(vel_basis_functions_y),
        s_vel_basis_function_z(vel_basis_functions_z), s_weight(weight),
        s_base_indices(base_indices) {}

  template <typename T> bool operator()(const T *coeffs_ak, T *residual) const {
    // Map coefficients to interpret the memory as eigen vector without copying
    Eigen::Map<const Eigen::VectorX<T>> coeffs(coeffs_ak,
                                               MAX_NUMBER_OF_VELOCITY_BASIS);

    // Create new vectors of converted types (templated types)
    Eigen::VectorX<T> t_vel_basis_functions_x =
        s_vel_basis_function_x.cast<T>();
    Eigen::VectorX<T> t_vel_basis_functions_y =
        s_vel_basis_function_y.cast<T>();
    Eigen::VectorX<T> t_vel_basis_functions_z =
        s_vel_basis_function_z.cast<T>();

    Eigen::Vector3<T> src_pt = utilities::toTemplatedPoint<T>(s_source_point);
    Eigen::Vector3<T> target_pt =
        utilities::toTemplatedPoint<T>(s_target_point);

    Eigen::Vector3<T> updated_pos = src_pt;
    PositionIncrementor<T> position_incrementor(&updated_pos);
    position_incrementor.incrementPosition(
        t_vel_basis_functions_x, t_vel_basis_functions_y,
        t_vel_basis_functions_z, coeffs, s_base_indices, NUMBER_OF_TIME_STEPS);

    residual[0] = (target_pt - updated_pos).norm() * T(s_weight);

    return true;
  }

  const Eigen::Vector3d &s_source_point;
  const Eigen::Vector3d &s_target_point;
  const Eigen::VectorXd &s_vel_basis_function_x;
  const Eigen::VectorXd &s_vel_basis_function_y;
  const Eigen::VectorXd &s_vel_basis_function_z;
  const double s_weight;
  const std::vector<adi::deformation_field::BasisIndices> &s_base_indices;
};

struct GaussianPriorCostFunctor {
  GaussianPriorCostFunctor(const Eigen::MatrixXd &L_inv) : s_L_inv(L_inv) {}

  template <typename T> bool operator()(const T *coeffs_ak, T *residual) const {
    // Compute the Gaussian prior residual term with (1/2) * (a^T * L_inv
    //* a)
    Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> coeffs(
        coeffs_ak, s_L_inv.rows());
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> s_L_inv_T =
        s_L_inv.template cast<T>();

    T quadratic_term = coeffs.transpose() * s_L_inv_T * coeffs;
    residual[0] = T(0.5) * quadratic_term;

    return true;
  }

  const Eigen::MatrixXd &s_L_inv; // Reference to the L_inv matrix
};

void Optimize(
    const Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
                         MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
        &vel_basis_functions,
    Eigen::VectorXd *coeffs_ak, const std::vector<adi::Point> &target_cloud,
    const std::vector<adi::Point> &src_cloud,
    const Eigen::MatrixXd &soft_corr_matrix, const Eigen::MatrixXd &L_inv,
    const std::vector<adi::deformation_field::BasisIndices> &base_indices) {
  // Prepare the ceres problem
  ceres::Problem problem;

  for (size_t i = 0; i < src_cloud.size(); ++i) {
    Eigen::Vector3d src_point = src_cloud.at(i).s_point;
    adi::deformation_field::DeformationField df;
    auto [vel_basis_function_x, vel_basis_function_y, vel_basis_function_z] =
        df.computeVelocityBasisPerPointPerDimension<double>(
            vel_basis_functions(i, Fastor::all, Fastor::all));
    for (size_t j = 0; j < target_cloud.size(); ++j) {
      Eigen::Vector3d target_point = target_cloud.at(j).s_point;
      const double weight = soft_corr_matrix(i, j);
      problem.AddResidualBlock(
          new ceres::AutoDiffCostFunction<MatchingCostFunctor, 1,
                                          MAX_NUMBER_OF_VELOCITY_BASIS>(
              new MatchingCostFunctor(
                  src_point, target_point, vel_basis_function_x,
                  vel_basis_function_y, vel_basis_function_z, weight,
                  base_indices)),
          nullptr, coeffs_ak->data());
    }
  }

  // Add Gaussian prior residual block
  problem.AddResidualBlock(
      new ceres::AutoDiffCostFunction<GaussianPriorCostFunctor, 1,
                                      MAX_NUMBER_OF_VELOCITY_BASIS>(
          new GaussianPriorCostFunctor(L_inv)),
      nullptr,
      coeffs_ak->data() // Optimizing coeffs_ak
  );

  // Set solver options
  ceres::Solver::Options options;
  options.max_num_iterations = 10;
  options.linear_solver_type = ceres::DENSE_QR;

  // Solve the problem
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  // Print solver summary
  std::cout << summary.FullReport() << std::endl;
}

} // namespace numerics
} // namespace adi
