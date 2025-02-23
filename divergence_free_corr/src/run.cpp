#include <chrono>

#include "../external/include/clue.hpp"
#include "parallel_optimizer.hpp"
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
           << source_downsampled_cloud.size());

  LOG_INFO("Size of target cloud for processing: "
           << target_downsampled_cloud.size());

  // Compute initial correspondence based on SHOT features
  std::unique_ptr<std::vector<std::pair<adi::Point, adi::Point>>>
      initial_correspondences = utilities::computeCorrespondences(
          source_downsampled_cloud, target_downsampled_cloud);
  LOG_INFO("Initial correspondences computed");

  uint32_t iter = 0;
  std::array<adi::deformation_field::BasisIndices, MAX_NUMBER_OF_VELOCITY_BASIS>
      base_indices =
          adi::numerics::GenerateBasisIndices(MAX_NUMBER_OF_VELOCITY_BASIS);
  LOG_INFO("Basis indices generated ");

  // Compute covariance matrix of gaussian distribution
  Eigen::MatrixXd L_inv = adi::numerics::computePrecisionMatrix(base_indices);
  LOG_INFO("Completed computing precision matrix");

  auto start = std::chrono::high_resolution_clock::now();
  // Update point cloud Y with the current deformation field
  while (iter < MAX_NUMBER_OF_ITERS) {
    // E-step : Compute soft correspondences
    auto matches = adi::matching::Matching(*initial_correspondences);
    Eigen::MatrixXd soft_corr_matrix = matches.computeSoftCorrespondences(
        source_downsampled_cloud, target_downsampled_cloud);

    LOG_INFO("Soft correspondences computed");

    // M-step : handled by ceres

    adi::parallel_optimizer::ParallelOptimizer optimizer(MAX_NUMBER_OF_THREADS);
    optimizer.Optimize(target_downsampled_cloud, source_downsampled_cloud,
                       soft_corr_matrix, L_inv, base_indices);
    iter = iter + 1;
    // updated_cloud.clear();
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  {
    std::stringstream msg;
    msg << "Duration of optimization (in seconds): " << duration.count()
        << std::endl;
    LOG_INFO(msg.str());
  }
  return;
}