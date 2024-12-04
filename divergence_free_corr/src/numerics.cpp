#include "numerics.hpp"
#include "../external/include/clue.hpp"
#include "utilities.hpp"

#include <Fastor/Fastor.h>
#include <ceres/autodiff_cost_function.h>
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
  Eigen::MatrixXd precision_matrix = Eigen::MatrixXd::Identity(
      MAX_NUMBER_OF_VELOCITY_BASIS, MAX_NUMBER_OF_VELOCITY_BASIS);

  Eigen::VectorXd inverse_eigen_vals =
      Eigen::VectorXd::Zero(MAX_NUMBER_OF_VELOCITY_BASIS);

  for (uint32_t i = 0; i < MAX_NUMBER_OF_VELOCITY_BASIS; ++i) {
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
void PositionIncrementor<T>::computeVelocityField(
    const Eigen::VectorXd &vel_basis_x, const Eigen::VectorXd &vel_basis_y,
    const Eigen::VectorXd &vel_basis_z, T *vel_field) {

  T dot_x = this->dot(vel_basis_x);
  T dot_y = this->dot(vel_basis_y);
  T dot_z = this->dot(vel_basis_z);

  vel_field[0] = dot_x;
  vel_field[1] = dot_y;
  vel_field[2] = dot_z;

  return;
}

template <typename T>
void PositionIncrementor<T>::incrementPosition(
    T *to_be_updated_pt, const Eigen::VectorXd &vel_basis_x,
    const Eigen::VectorXd &vel_basis_y, const Eigen::VectorXd &vel_basis_z,
    const std::vector<adi::deformation_field::BasisIndices> &basis_indices,
    const double num_time_steps) {

  const T dt = T(1.0) / T(num_time_steps);

  for (uint32_t t = 0; t < num_time_steps; ++t) {
    T intermediate_vel_field[3];
    utilities::fillVector(intermediate_vel_field, Eigen::Vector3d::Zero());
    this->computeVelocityField(vel_basis_x, vel_basis_y, vel_basis_z,
                               intermediate_vel_field);
    std::cout << "intermediate_vel_field[0]: " << intermediate_vel_field[0]
              << std::endl;
    to_be_updated_pt[0] += T(0.5) * dt * intermediate_vel_field[0];
    std::cout << "updated_pt[0]: " << to_be_updated_pt[0] << std::endl;
    to_be_updated_pt[1] += T(0.5) * dt * intermediate_vel_field[1];
    to_be_updated_pt[2] += T(0.5) * dt * intermediate_vel_field[2];

    // adi::deformation_field::DeformationField df;
    // Eigen::Matrix<T, MAX_NUMBER_OF_VELOCITY_BASIS, 3>
    //     vel_basis_functions_updated =
    //         df.computeVelocityBasisFunctionsPerPoint<T>(basis_indices,
    //                                                     *s_updated_pt);

    // Eigen::Vector3<T> final_vel_field = this->computeVelocityField(
    //     vel_basis_functions_updated.col(0),
    //     vel_basis_functions_updated.col(1),
    //     vel_basis_functions_updated.col(2), coeffs);

    // (*s_updated_pt)[0] += dt * final_vel_field[0];
    // (*s_updated_pt)[1] += dt * final_vel_field[1];
    // (*s_updated_pt)[2] += dt * final_vel_field[2];
  }
}

struct MatchingCostFunctor {
  MatchingCostFunctor(
      const Eigen::Vector3d &source_point, const Eigen::Vector3d &target_point,
      const double weight,
      const std::vector<adi::deformation_field::BasisIndices> &base_indices)
      : s_target_point(target_point), s_source_point(source_point),
        s_weight(weight), s_base_indices(base_indices) {

    // TODO: compute velocity basis function here for each point
    adi::deformation_field::DeformationField s_df;
    Eigen::MatrixXd vel_basis_functions =
        s_df.computeVelocityBasisFunctionsPerPoint<double>(base_indices,
                                                           source_point);

    s_vel_basis_function_x = vel_basis_functions.col(0);
    s_vel_basis_function_y = vel_basis_functions.col(1);
    s_vel_basis_function_z = vel_basis_functions.col(2);

    if (s_vel_basis_function_x.hasNaN() || s_vel_basis_function_y.hasNaN() ||
        s_vel_basis_function_z.hasNaN()) {
      std::cout << "Nan values present in velocity basis functions"
                << std::endl;
    }
  }

  template <typename T> bool checkForNaN(const Eigen::VectorX<T> &vec) const {
    for (int i = 0; i < vec.size(); ++i) {
      if (ceres::isnan(vec[i])) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool operator()(const T *const coeffs_ak, T *residual) const {
    T src_pt[3];
    T target_pt[3];
    utilities::toTemplatedPoint(s_source_point, src_pt);
    utilities::toTemplatedPoint(s_target_point, target_pt);

    PositionIncrementor<T> position_incrementor((T *const)(coeffs_ak));
    position_incrementor.incrementPosition(
        src_pt, s_vel_basis_function_x, s_vel_basis_function_y,
        s_vel_basis_function_z, s_base_indices, NUMBER_OF_TIME_STEPS);

    // if (checkForNaN<T>(updated_pos)) {
    //   std::cout << "Updated position is invalid" << std::endl;
    //   return false;
    // }

    // residual[0] = (s_target_point - updated_pos).norm() * T(s_weight);
    residual[0] = T(s_weight);
    // std::cout << "Face value of the residual: " << residual[0] << std::endl;
    return true;
  }

  const Eigen::Vector3d &s_source_point;
  const Eigen::Vector3d &s_target_point;
  Eigen::VectorXd s_vel_basis_function_x;
  Eigen::VectorXd s_vel_basis_function_y;
  Eigen::VectorXd s_vel_basis_function_z;
  const double s_weight;
  const std::vector<adi::deformation_field::BasisIndices> &s_base_indices;
};

// struct GaussianPriorCostFunctor {
//   GaussianPriorCostFunctor(const Eigen::MatrixXd &L_inv) : s_L_inv(L_inv) {}

//   template <typename T> T computeQuadraticTerm(const T *const v2) const {
//     T sum = T(0.0);
//     for (uint32_t i = 0; i < MAX_NUMBER_OF_VELOCITY_BASIS; ++i) {
//       sum += T(s_L_inv(i, i)) * v2[i] * v2[i];
//     }

//     return sum;
//   }

//   template <typename T> bool operator()(const T *coeffs_ak, T *residual)
//   const {
//     // Compute the Gaussian prior residual term with (1/2) * (a^T * L_inv
//     //* a)
//     Eigen::Map<const Eigen::Matrix<T, MAX_NUMBER_OF_VELOCITY_BASIS, 1>>
//     coeffs(
//         coeffs_ak, MAX_NUMBER_OF_VELOCITY_BASIS);

//     T quadratic_term = this->computeQuadraticTerm<T>(coeffs.data());
//     residual[0] = T(0.5) * quadratic_term;

//     return true;
//   }

//   const Eigen::MatrixXd &s_L_inv; // Reference to the L_inv matrix
// };

void Optimize(
    const std::vector<adi::Point> &target_cloud,
    const std::vector<adi::Point> &src_cloud,
    const Eigen::MatrixXd &soft_corr_matrix, const Eigen::MatrixXd &L_inv,
    const std::vector<adi::deformation_field::BasisIndices> &base_indices) {
  // Prepare the ceres problem
  ceres::Problem problem;

  double coeffs_ak[MAX_NUMBER_OF_VELOCITY_BASIS];
  utilities::fillVector<double>(
      coeffs_ak, Eigen::VectorXd::Zero(MAX_NUMBER_OF_VELOCITY_BASIS));
  PositionIncrementor<double> position_incrementor(coeffs_ak);

  for (size_t i = 0; i < src_cloud.size(); ++i) {
    Eigen::Vector3d src_point = src_cloud.at(i).s_point;
    for (size_t j = 0; j < target_cloud.size(); ++j) {
      Eigen::Vector3d target_point = target_cloud.at(j).s_point;

      // TODO: WEIGHTS CAN BE APPLIED TO THE LOSS
      const double weight = soft_corr_matrix(i, j);
      problem.AddResidualBlock(
          new ceres::AutoDiffCostFunction<MatchingCostFunctor, 1,
                                          MAX_NUMBER_OF_VELOCITY_BASIS>(
              new MatchingCostFunctor(src_point, target_point, weight,
                                      base_indices)),
          nullptr, position_incrementor.getCoeffs());
    }
  }

  // Add Gaussian prior residual block
  // problem.AddResidualBlock(
  //     new ceres::AutoDiffCostFunction<GaussianPriorCostFunctor, 1,
  //                                     MAX_NUMBER_OF_VELOCITY_BASIS>(
  //         new GaussianPriorCostFunctor(L_inv)),
  //     nullptr,
  //     coeffs_ak->data() // Optimizing coeffs_ak
  // );

  // Set solver options
  ceres::Solver::Options options;
  options.max_num_iterations = 100;
  options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;

  // Solve the problem
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  // Print solver summary
  std::cout << summary.FullReport() << std::endl;
}

} // namespace numerics
} // namespace adi
