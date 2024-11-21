#include "numerics.hpp"
#include "../external/include/clue.hpp"
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

// template <typename T>
// Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3> RungeKutta2Integration(
// const std::vector<adi::deformation_field::BasisIndices>& basis_indices,
// const Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3>& src_point_cloud,
// const Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 1>& coeffs_ak,
// const Fastor::Tensor<T, NUMBER_OF_SAMPLE_POINTS,
// MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>& vel_basis_functions, const
// uint32_t num_time_steps)
// {
//     Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3> updated_pts =
//     src_point_cloud;

//     const T dt = T(1.0) / T(num_time_steps);

//     for (uint32_t t = 0; t < num_time_steps; ++t) {
//         adi::deformation_field::DeformationField df;

//         // Compute initial velocity field at current point positions
//         Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3> intermediate_vel_field =
//             adi::deformation_field::DeformationField::computeVelocityFieldTemplated(coeffs_ak,
//             vel_basis_functions);

//         // Compute midpoint based on current velocity field
//         Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3> midpoint = updated_pts;
//         midpoint.col(0) += T(0.5) * dt * intermediate_vel_field.col(0);
//         midpoint.col(1) += T(0.5) * dt * intermediate_vel_field.col(1);
//         midpoint.col(2) += T(0.5) * dt * intermediate_vel_field.col(2);

//         // Update velocity basis functions based on the midpoint positions
//         auto updated_vel_basis_functions =
//         df.computeVelocityBasisFunctions(basis_indices, midpoint);

//         // Compute updated velocity field at midpoint positions
//         Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3> updated_vel_field =
//             adi::deformation_field::DeformationField::computeVelocityFieldTemplated(coeffs_ak,
//             updated_vel_basis_functions);

//         // Update points with final velocity field
//         updated_pts.col(0) += dt * updated_vel_field.col(0);
//         updated_pts.col(1) += dt * updated_vel_field.col(1);
//         updated_pts.col(2) += dt * updated_vel_field.col(2);
//     }

//     return updated_pts;
// }

const double EucledianDistance(const Eigen::Vector3d &point_1,
                               const Eigen::Vector3d &point_2) {
  return (point_1 - point_2).norm();
}

// struct MatchingCostFunctor{
//     MatchingCostFunctor(const Eigen::MatrixXd &target_point,
//                         const Eigen::MatrixXd &src_cloud,
//                         const Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
//                         MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
//                         &vel_basis_functions, const Eigen::VectorXd weight,
//                         const
//                         std::vector<adi::deformation_field::BasisIndices>
//                         &base_indices) :s_target_point(target_point),
//                          s_src_cloud(src_cloud),
//                          s_vel_basis_functions(vel_basis_functions),
//                          s_weight(weight),
//                          s_base_indices(base_indices){}

//     template <typename T>
//     bool operator()(const T* coeffs_ak, T* residual) const{
//         // Map coefficients
//         Eigen::Map<const Eigen::Matrix<T, MAX_NUMBER_OF_VELOCITY_BASIS, 1>>
//         coeffs(coeffs_ak, MAX_NUMBER_OF_VELOCITY_BASIS);

//         // Map target cloud points to templated type
//         Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3> target_point =
//         s_target_point.template cast<T>();

//         // Map source points to templated type
//         Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3> src_points =
//         s_src_cloud.template cast<T>();

//         Fastor::Tensor<T, NUMBER_OF_SAMPLE_POINTS,
//         MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
//         templated_vel_basis_functions(s_vel_basis_functions);

//         auto new_src_cloud =
//         RungeKutta2IntegrationTemplated<T>(s_base_indices, src_points,
//         coeffs, templated_vel_basis_functions, 10);

//         Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 1> residuals =
//         utilities::L2Norm<T>(new_src_cloud, target_point);

//         // Compute weighted L2 norms
//         Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 1> weighted_residuals =
//         residuals * s_weight.template cast<T>().array();

//         // Sum the weighted residuals
//         T total_weighted_residual = weighted_residuals.sum();

//         residual[0] = total_weighted_residual;

//         return true;
//     }

//     const Eigen::MatrixXd &s_target_point;
//     const Eigen::MatrixXd &s_src_cloud;
//     const Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
//     MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH> s_vel_basis_functions; const
//     Eigen::VectorXd s_weight; const
//     std::vector<adi::deformation_field::BasisIndices> &s_base_indices;
// };

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

  // for (size_t i = 0; i < target_cloud.size(); ++i) {
  //   // boradcast the target cloud
  //   Eigen::MatrixXd broadcasted_target_point =
  //       Eigen::MatrixXd::Ones(NUMBER_OF_SAMPLE_POINTS, 3);
  //   broadcasted_target_point.col(0).setConstant(target_cloud.at(i).s_point.x());
  //   broadcasted_target_point.col(1).setConstant(target_cloud.at(i).s_point.y());
  //   broadcasted_target_point.col(2).setConstant(target_cloud.at(i).s_point.z());

  //   auto weight = soft_corr_matrix.col(i);

  //   // Create a residual block for this target point
  //   problem.AddResidualBlock(
  //       new ceres::AutoDiffCostFunction<MatchingCostFunctor, 1,
  //                                       MAX_NUMBER_OF_VELOCITY_BASIS>(
  //           new MatchingCostFunctor(broadcasted_target_point,
  //                                   utilities::toEigenMatrix(src_cloud),
  //                                   vel_basis_functions, weight,
  //                                   base_indices)),
  //       nullptr, coeffs_ak->data());
  // }

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
