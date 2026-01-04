#ifndef FEATURES_HPP
#define FEATURES_HPP

#include <pcl/features/shot_omp.h>

namespace adi {
namespace features {
/**
 * @brief Function to compute SHOT features
 *
 * @param cloud_with_normals
 * @param search_radius
 *
 * @return Pointer to the computed SHOT features
 */
pcl::PointCloud<pcl::SHOT352>::Ptr computeShotFeatures(
    const pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals,
    const double search_radius);
} // namespace features
} // namespace adi
#endif