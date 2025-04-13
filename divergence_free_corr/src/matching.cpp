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

      if (isnan(metric_distance(i, j)))
        metric_distance(i, j) = 0.0;
    }
  }

  return metric_distance;
}

Eigen::SparseMatrix<double> Matching::computeSoftCorrespondences(
    const std::vector<adi::Point> &source_point_cloud,
    const std::vector<adi::Point> &target_point_cloud) {

  const double scaling_factor = this->ComputeScalingFactor();

  const Eigen::MatrixXd metric_distance = computeMetricDistance(
      source_point_cloud, target_point_cloud, scaling_factor);

  // Create triplet list for sparse matrix construction
  std::vector<Eigen::Triplet<double>> triplets;
  // Reserve estimated capacity (adjust based on expected sparsity)
  triplets.reserve(source_point_cloud.size() * target_point_cloud.size() / 10);

  for (uint32_t i = 0; i < source_point_cloud.size(); ++i) {
    const double denominator = std::exp(metric_distance.row(i).array().sum() *
                                        (-1 / (2 * SIGMA * SIGMA))) +
                               std::pow((2 * M_PI * SIGMA * SIGMA), 1.5);
    for (uint32_t j = 0; j < target_point_cloud.size(); ++j) {
      double distance = metric_distance(i, j);
      const double numerator =
          std::exp((-1 / (2 * SIGMA * SIGMA) * distance * distance));

      double value = numerator / denominator;

      if (!std::isnan(value) && std::abs(value) > 1e-10)
        triplets.emplace_back(i, j, value);
    }
  }

  // Create and populate the sparse matrix
  Eigen::SparseMatrix<double> correspondences(source_point_cloud.size(),
                                              target_point_cloud.size());
  correspondences.setFromTriplets(triplets.begin(), triplets.end());

  // Compress the sparse matrix
  correspondences.makeCompressed();

  return correspondences;
}
} // namespace matching
} // namespace adi