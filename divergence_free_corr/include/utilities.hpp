#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "deformation_field.hpp"
#include "io.hpp"

#include <ceres/jet.h>
#include <random>

namespace utilities {
/**
 * @brief Compute eucledian distance between 2 points
 *
 * @param point_1
 * @param point_2
 *
 * @return Eucledian distance between 2 points
 */
double eucledianDistance(const Eigen::Vector3d &point_1,
                         const Eigen::Vector3d &point_2);

/**
 * @brief Computes L2 norm between 2 descriptors
 *
 * @param descriptor_1
 * @param descriptor_2
 *
 * @return L2 norm between 2 descriptors
 */
double computeL2Norm(const std::array<double, 352> &descriptor_1,
                     const std::array<double, 352> &descriptor_2);

/**
 * @brief Checks if the feature computed is valid
 *
 * @param descriptor
 *
 * @return A boolean if the feature computed is valid
 */
bool isValidDescriptor(const pcl::SHOT352 &descriptor);

/**
 * @brief Computes the initial set of hard correspondences
 *
 * @param source_point_cloud
 * @param target_point_cloud
 *
 * @return Pointer to the hard correspondences computed
 */
std::unique_ptr<std::vector<std::pair<adi::Point, adi::Point>>>
computeCorrespondences(std::vector<adi::Point> source_point_cloud,
                       std::vector<adi::Point> target_point_cloud);

/**
 * @brief Converts to eigen matrix
 *
 * @param cloud
 *
 * @return Point cloud wherein each column represents values in each dimension
 */
Eigen::MatrixXd toEigenMatrix(const std::vector<adi::Point> &cloud);

/**
 * @brief Convert to vector of eigen 3D vectors
 *
 * @param cloud
 * @param point_cloud
 */
std::vector<adi::Point> toPointCloud(const Eigen::MatrixXd &cloud);

const Eigen::MatrixXd computeLInv(
    const std::vector<adi::deformation_field::BasisIndices> &basis_indices);

template <typename T>
inline void toTemplatedPoint(const Eigen::Vector3d &pt, T *output) {
  output[0] = T(pt.x());
  output[1] = T(pt.y());
  output[2] = T(pt.z());

  return;
}

template <typename T> void fillVector(T *v1, const Eigen::VectorXd &v2) {
  for (uint32_t i = 0; i < v2.size(); ++i) {
    v1[i] = T(v2[i]);
  }
  return;
}

template <typename T>
void fillArray(T *v1,
               const std::array<adi::deformation_field::BasisIndices,
                                MAX_NUMBER_OF_VELOCITY_BASIS> &basis_indices) {
  // random number generator
  std::random_device rd;
  std::mt19937 gen(rd());

  for (uint32_t i = 0; i < MAX_NUMBER_OF_VELOCITY_BASIS; ++i) {
    std::normal_distribution<double> dist(
        0.0, std::sqrt(basis_indices[i].eigen_value));
    v1[i] = dist(gen);
  }
}

} // namespace utilities
#endif // UTILITIES_HPP