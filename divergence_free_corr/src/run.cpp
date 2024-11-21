#include <chrono>

#include "../external/include/clue.hpp"
#include "run.hpp"

void run(adi::pointCloud *source_cloud, adi::pointCloud *target_cloud) {
  // downsample the point clouds
  source_cloud->samplePointCloud(NUMBER_OF_SAMPLE_POINTS);
  target_cloud->samplePointCloud(NUMBER_OF_SAMPLE_POINTS);

  // Compute features
  source_cloud->computeShotFeatures(SEARCH_RADIUS);
  LOG_INFO("Completed computing SHOT features for source cloud");

  target_cloud->computeShotFeatures(SEARCH_RADIUS);
  LOG_INFO("Completed computing SHOT features for target cloud");

  // get downsampled cloud with features
  std::vector<adi::Point> source_downsampled_cloud =
      source_cloud->getPointCloud();
  std::vector<adi::Point> target_downsampled_cloud =
      target_cloud->getPointCloud();

  LOG_INFO("Size of source cloud for processing: "
           << source_downsampled_cloud.size() << "\n"
           << "Size of target cloud for processing: "
           << target_downsampled_cloud.size());

  // Serialize downsampled source cloud
  // adi::pointCloud::serializeCloud(source_downsampled_cloud,
  // "/workspaces/divergence_free_correspondence/data/source_downsampled_cloud.ply");
  // std::cout << "Completed writing source downsampled cloud" << std::endl;

  // Compute initial correspondence based on SHOT features
  std::unique_ptr<std::vector<std::pair<adi::Point, adi::Point>>>
      initial_correspondences = utilities::computeCorrespondences(
          source_downsampled_cloud, target_downsampled_cloud);
  LOG_INFO("Initial correspondences computed");

  // Gauss Newton optimization looop until convergence
  uint32_t iter = 0;
  std::vector<adi::deformation_field::BasisIndices> base_indices =
      adi::numerics::GenerateBasisIndices(MAX_NUMBER_OF_VELOCITY_BASIS);
  LOG_INFO("Basis indices generated ");

  // precompute the velocity basis functions as a matrix
  adi::deformation_field::DeformationField df;
  auto vel_basis_functions = df.computeVelocityBasisFunctions(
      base_indices, utilities::toEigenMatrix(source_downsampled_cloud));

  LOG_INFO("Completed computing velocity basis functions ");
  Eigen::VectorXd coeffs_ak =
      Eigen::VectorXd::Zero(MAX_NUMBER_OF_VELOCITY_BASIS);

  // Update point cloud Y with the current deformation field
  while (iter < MAX_NUMBER_OF_ITERS) {
    // RK2 integration
    auto start = std::chrono::high_resolution_clock::now();
    Eigen::MatrixXd updated_pts = adi::numerics::RungeKutta2Integration(
        base_indices, utilities::toEigenMatrix(source_downsampled_cloud),
        coeffs_ak, vel_basis_functions, NUMBER_OF_TIME_STEPS);
    auto end = std::chrono::high_resolution_clock::now();

    std::vector<adi::Point> updated_cloud =
        utilities::toPointCloud(updated_pts);

    std::chrono::duration<double, std::milli> duration = end - start;
    LOG_INFO("Time taken (in milliseconds) for "
             << NUMBER_OF_TIME_STEPS
             << " time steps of RK2 integration: " << duration.count());

    LOG_INFO(" Size of updated cloud: " << updated_cloud.size());
    LOG_INFO(" Completed RK integration for iteration: " << iter);

    // Serialize downsampled source cloud
    // adi::pointCloud::serializeCloud(updated_cloud,
    // "/workspaces/divergence_free_correspondence/data/source_cloud_updated.ply");
    // std::cout << "Completed writing updated source cloud" << std::endl;

    // E-step : Compute soft correspondences
    auto matches = adi::matching::Matching(*initial_correspondences);
    Eigen::MatrixXd soft_corr_matrix = matches.computeSoftCorrespondences(
        updated_cloud, target_downsampled_cloud);

    LOG_INFO("Soft correspondences computed");
    // Compute updated velocity basis functions of the updated cloud
    // auto vel_basis_functions_updated =
    //     df.computeVelocityBasisFunctions(base_indices, updated_pts);
    // LOG_INFO("Completed computing updated velocity field");
    // M-step : handled by ceres

    //     // M-step : Gauss Newton optimization
    //     Eigen::MatrixXd hessian_matrix =
    //     adi::numerics::ComputeHessian(base_indices, updated_cloud,
    //     target_downsampled_cloud, soft_corr_matrix, r0);

    //     Eigen::VectorXd gradient =
    //     adi::numerics::ComputeGradient(base_indices, updated_cloud,
    //     target_downsampled_cloud, soft_corr_matrix, r0);

    //     Eigen::VectorXd delta_a = -hessian_matrix.ldlt().solve(gradient);

    //     coeffs_ak += delta_a;

    //     std::cout << "Updated coefficients of deformation field for iteration
    //     " << iter << std::endl; if(delta_a.norm() < TOLERANCE)
    //     {
    //         std::cout << "Coefficients of deformation field converged" <<
    //         std::endl; break;
    //     }
    iter = iter + 1;
    updated_cloud.clear();
  }

  // for(size_t a_k_i = 0;a_k_i < coeffs_ak.size();++a_k_i)
  // {
  //     std::cout << coeffs_ak[a_k_i] << std::endl;
  // }
  return;
}