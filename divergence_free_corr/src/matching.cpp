#include "matching.hpp"
#include "constants.hpp"
#include "utilities.hpp"

namespace adi {
namespace matching {

Matching::Matching(
    const std::vector<std::pair<adi::Point, adi::Point>> &correspondences)
    : m_correspondences(correspondences) {}

const double Matching::computeMeanEucledianDistance() {
  double mean_eucledian_distance = 0;
  for (uint32_t i = 0; i < m_correspondences.size(); ++i) {
    mean_eucledian_distance += (m_correspondences[i].first.s_point -
                                m_correspondences[i].second.s_point)
                                   .norm();
  }
  return mean_eucledian_distance / m_correspondences.size();
}

const double Matching::computeMeanDescriptorDistance() {
  double mean_descriptor_distance = 0;
  for (uint32_t i = 0; i < m_correspondences.size(); ++i) {
    mean_descriptor_distance +=
        utilities::computeL2Norm(m_correspondences[i].first.s_descriptor,
                                 m_correspondences[i].second.s_descriptor);
  }
  return mean_descriptor_distance / m_correspondences.size();
}

const double Matching::ComputeScalingFactor() {
  const double mean_eucledian_distance = this->computeMeanEucledianDistance();
  const double mean_descriptor_distance = this->computeMeanDescriptorDistance();
  const double factor = mean_eucledian_distance / mean_descriptor_distance;

  return factor;
}

const Eigen::MatrixXd Matching::computeMetricDistance(
    const std::vector<adi::Point> &source_point_cloud,
    const std::vector<adi::Point> &target_point_cloud,
    const double &scaling_factor) {
  Eigen::MatrixXd metric_distance = Eigen::MatrixXd::Zero(
      source_point_cloud.size(), source_point_cloud.size());

  for (uint32_t i = 0; i < source_point_cloud.size(); ++i) {
    for (uint32_t j = 0; j < target_point_cloud.size(); ++j) {
      const double local_eucledian_distance =
          (source_point_cloud[i].s_point - target_point_cloud[j].s_point)
              .norm();
      const double local_descriptor_distance =
          utilities::computeL2Norm(source_point_cloud[i].s_descriptor,
                                   target_point_cloud[j].s_descriptor);
      metric_distance(i, j) =
          local_eucledian_distance + scaling_factor * local_descriptor_distance;

      std::cout << "Metric distance between points: " << i << " and " << j
                << " : " << metric_distance(i, j) << std::endl;
    }
  }

  return metric_distance;
}

Eigen::MatrixXd Matching::computeSoftCorrespondences(
    const std::vector<adi::Point> &source_point_cloud,
    const std::vector<adi::Point> &target_point_cloud) {
  Eigen::MatrixXd correspondences = Eigen::MatrixXd::Zero(
      source_point_cloud.size(), source_point_cloud.size());
  const double scaling_factor = this->ComputeScalingFactor();
  std::cout << "Scaling factor: " << scaling_factor << std::endl;
  const Eigen::MatrixXd metric_distance = computeMetricDistance(
      source_point_cloud, target_point_cloud, scaling_factor);

  for (uint32_t i = 0; i < source_point_cloud.size(); ++i) {
    for (uint32_t j = 0; j < target_point_cloud.size(); ++j) {
      double distance = metric_distance(i, j);

      const double numerator =
          std::exp((-1 / (2 * SIGMA * SIGMA) * distance * distance));
      const double denominator = std::exp(metric_distance.row(i).array().sum() *
                                          (-1 / (2 * SIGMA * SIGMA))) +
                                 std::pow((2 * M_PI * SIGMA * SIGMA), 1.5);
      correspondences(i, j) = numerator / denominator;
    }
  }

  return correspondences;
}

// Eigen::MatrixXd Matching::computeSoftCorrespondences(
//     const std::vector<adi::Point> &source_point_cloud,
//     const std::vector<adi::Point> &target_point_cloud) {

//   Eigen::MatrixXd correspondences = Eigen::MatrixXd::Zero(
//       source_point_cloud.size(), source_point_cloud.size());
//   const Eigen::MatrixXd metric_distance =
//       computeMetricDistance(source_point_cloud, target_point_cloud);

//   const double sigma_squared = SIGMA * SIGMA;

//   // Debug print distances
//   // std::cout << "Distances range: " << metric_distance.minCoeff() << " to "
//   //           << metric_distance.maxCoeff() << std::endl;
//   // std::cout << "Sigma squared: " << sigma_squared << std::endl;

//   for (uint32_t i = 0; i < source_point_cloud.size(); ++i) {
//     // Compute all exp(-d²/2σ²) terms for this row
//     Eigen::VectorXd exp_terms(target_point_cloud.size());
//     for (uint32_t j = 0; j < target_point_cloud.size(); ++j) {
//       const double distance = metric_distance(i, j);
//       exp_terms(j) = std::exp(-distance * distance / (2 * sigma_squared));
//     }

//     // Compute denominator (sum of all exp terms)
//     double denominator = exp_terms.sum();

//     // Debug check
//     // std::cout << "Row " << i << " denominator: " << denominator <<
//     std::endl;

//     // Safety check to avoid division by zero
//     if (denominator < 1e-10) {
//       denominator = 1e-10;
//     }

//     // Compute final correspondences for this row
//     for (uint32_t j = 0; j < target_point_cloud.size(); ++j) {
//       correspondences(i, j) = exp_terms(j) / denominator;

//       // Debug check
//       // if (std::isnan(correspondences(i, j))) {
//       //   std::cout << "NaN at (" << i << "," << j
//       //             << "), exp_term: " << exp_terms(j)
//       //             << ", denominator: " << denominator << std::endl;
//       // }
//     }

//     // Debug check: row sum should be close to 1
//     // std::cout << "Row " << i << " sum: " << correspondences.row(i).sum()
//     //           << std::endl;
//   }

//   return correspondences;
// }
} // namespace matching
} // namespace adi