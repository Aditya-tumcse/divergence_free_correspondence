#ifndef MATCHING_HPP
#define MATCHING_HPP

#include "io.hpp"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SparseCore>

namespace adi {
namespace matching {

class Matching {
public:
  Matching(
      const std::vector<std::pair<adi::Point, adi::Point>> &correspondences);

  /**
   * @brief Compute soft correspondences that follow GMM
   *
   * @param source_point_cloud
   * @param target_point_cloud
   *
   * @return Eigen matrix of soft correspondences between source and target
   * pointcloud
   */
  Eigen::SparseMatrix<double>
  computeSoftCorrespondences(const std::vector<adi::Point> &source_point_cloud,
                             const std::vector<adi::Point> &target_point_cloud);

  ~Matching() = default;

private:
  const std::vector<std::pair<adi::Point, adi::Point>> &m_correspondences;

  const double computeMeanEucledianDistance();

  const double computeMeanDescriptorDistance();

  const double ComputeScalingFactor();

  const Eigen::MatrixXd
  computeMetricDistance(const std::vector<adi::Point> &source_point_cloud,
                        const std::vector<adi::Point> &target_point_cloud,
                        const double &scaling_factor);
};
} // namespace matching
} // namespace adi

#endif