#ifndef IO_HPP
#define IO_HPP

#include <Eigen/Dense>
#include <array>
#include <pcl/common/io.h>
#include <pcl/features/normal_3d.h>
#include <pcl/io/ply_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include "constants.hpp"

namespace adi {

/**
 * @brief Structure for storing information regarding the point of the
 * pointcloud
 *
 */
struct Point {
  Eigen::Vector3d s_point;              // Point Coordinates
  Eigen::Vector3d s_normal;             // Normal info
  std::array<double, 352> s_descriptor; // Feature descriptor

  Point() : s_point(0.0, 0.0, 0.0), s_normal(0.0, 0.0, 0.0), s_descriptor{} {}

  Point(Eigen::Vector3d point, Eigen::Vector3d normal,
        std::array<double, 352> descriptor)
      : s_point(point), s_normal(normal), s_descriptor(descriptor) {}

  Point(pcl::PointXYZ point) : s_normal(0.0, 0.0, 0.0), s_descriptor{} {
    s_point.x() = point.x;
    s_point.y() = point.y;
    s_point.z() = point.z;
  }

  /**
   * @brief Convert to PCL format of the descriptor
   *
   * @return SHOT descriptor pcl format
   */
  pcl::SHOT352 convertToPCLDescriptor();
};

class pointCloud {
public:
  pointCloud(const std::string &point_cloud_path);

  pointCloud(const pcl::PointCloud<pcl::PointXYZ>::Ptr &source_cloud);

  /**
   * @brief Load the pointcloud file from disk
   *
   * @param cloud_file_path
   *
   * @return Pointer to the loaded pointcloud
   */
  pcl::PointCloud<pcl::PointXYZ>::Ptr
  loadPlyFile(const std::string &cloud_file_path);

  /**
   * @brief Set the loaded pointcloud from the disk into the format of this
   * class
   *
   * @param cloud
   */
  void setPointCloud(const std::unique_ptr<std::vector<Point>> cloud);

  /**
   * @brief Compute features on the pointcloud
   *
   * @param search_radius
   */
  void computeShotFeatures(const double &search_radius);

  /**
   * @brief Get the private member variable m_point_cloud
   *
   * @return The private memeber variable m_point_cloud of the class
   */
  const std::vector<Point> getPointCloud() const { return m_point_cloud; }

  /**
   * @brief Downsample the pointcloud
   *
   * @param max_number_of_points
   */
  void samplePointCloud(const uint32_t max_number_of_points);

  /**
   * @brief Get the size of the pointcloud
   *
   * @return Size of the pointcloud
   */
  const uint32_t getSize() const { return m_point_cloud.size(); }

  static void serializeCloud(const std::vector<adi::Point> &cloud,
                             const std::string &file_path);

  ~pointCloud() = default;

private:
  std::vector<Point> m_point_cloud;

  /**
   * @brief Convert to PCL format pointcloud
   *
   * @return Returns pointer to the PCL format pointcloud
   */
  pcl::PointCloud<pcl::PointXYZ>::Ptr toPclCloud();

  static pcl::PointCloud<pcl::PointXYZ>::Ptr
  toPclCloud(const std::vector<adi::Point> &cloud);

  /**
   * @brief Loads the pointcloud from PCL format to the format of the class
   *
   * @param point_cloud
   *
   */
  void loadPointCloud(const pcl::PointCloud<pcl::PointXYZ>::Ptr &point_cloud);

  /**
   * @brief Estimates normals of the pointcloud
   *
   * @param point_cloud
   * @param search_radius
   *
   * @return Pointer to the pointcloud in PCL format with estimated normals
   */
  pcl::PointCloud<pcl::PointNormal>::Ptr
  normalEstimation(pcl::PointCloud<pcl::PointXYZ>::Ptr point_cloud,
                   const double &search_radius);
};
} // namespace adi
#endif // IO_HPP