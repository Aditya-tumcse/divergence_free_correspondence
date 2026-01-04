#include "io.hpp"
#include "features.hpp"

#include <cassert>
#include <cmath>
#include <pcl/features/normal_3d_omp.h>

namespace adi {

pcl::SHOT352 Point::convertToPCLDescriptor() {
  pcl::SHOT352 descriptor;
  for (uint32_t i = 0; i < s_descriptor.size(); ++i) {
    descriptor.descriptor[i] = static_cast<float>(s_descriptor[i]);
  }
  return descriptor;
}

pointCloud::pointCloud(const std::string &point_cloud_path) {
  pcl::PointCloud<pcl::PointXYZ>::Ptr point_cloud =
      this->loadPlyFile(point_cloud_path);
  this->loadPointCloud(point_cloud);
}

void pointCloud::setPointCloud(
    const std::unique_ptr<std::vector<Point>> cloud) {
  m_point_cloud.clear();
  m_point_cloud.reserve(cloud->size());
  for (size_t i = 0; i < cloud->size(); ++i) {
    Point pt;
    pt.s_point.x() = cloud->at(i).s_point.x();
    pt.s_point.y() = cloud->at(i).s_point.y();
    pt.s_point.z() = cloud->at(i).s_point.z();

    m_point_cloud.push_back(pt);
  }
}

pointCloud::pointCloud(const pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud) {
  this->loadPointCloud(cloud);
}

pcl::PointCloud<pcl::PointXYZ>::Ptr
pointCloud::loadPlyFile(const std::string &cloud_file_path) {
  if (cloud_file_path.empty())
    return nullptr;
  else {
    pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud(
        new pcl::PointCloud<pcl::PointXYZ>);
    if (pcl::io::loadPLYFile<pcl::PointXYZ>(cloud_file_path, *input_cloud) ==
        -1) {
      std::cerr << "Error loading point cloud file" << std::endl;
      return nullptr;
    } else {
      return input_cloud;
    }
  }
}

void pointCloud::computeShotFeatures(const double &search_radius) {
  // convert downsampled cloud into pcl cloud
  pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_cloud = this->toPclCloud();

  // Clear the member variable
  m_point_cloud.clear();
  m_point_cloud.reserve(pcl_cloud->size());

  // Perform normal estimation
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals =
      this->normalEstimation(pcl_cloud, search_radius);

  if (cloud_with_normals->empty()) {
    std::cerr << "Error: Cloud with normals is empty. Exiting." << std::endl;
    std::exit(EXIT_FAILURE);
  } else {
    // Compute SHOT descriptors
    pcl::PointCloud<pcl::SHOT352>::Ptr shot_descriptor =
        adi::features::computeShotFeatures(cloud_with_normals, search_radius);

    // Check if sizes match
    if (cloud_with_normals->size() != shot_descriptor->size()) {
      std::cerr << "Error: Sizes of cloud with normals and SHOT descriptors do "
                   "not match."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // Store features
    for (std::size_t i = 0; i < cloud_with_normals->points.size(); ++i) {
      Eigen::Vector3d point, normal;
      std::array<double, 352> desc{};
      point[0] = cloud_with_normals->points[i].x;
      point[1] = cloud_with_normals->points[i].y;
      point[2] = cloud_with_normals->points[i].z;

      normal[0] = cloud_with_normals->points[i].normal_x;
      normal[1] = cloud_with_normals->points[i].normal_y;
      normal[2] = cloud_with_normals->points[i].normal_z;

      // Copy descriptor values
      for (std::size_t j = 0; j < shot_descriptor->at(i).descriptorSize();
           ++j) {
        desc[j] = shot_descriptor->at(i).descriptor[j];
      }

      // Store the point, normal, and descriptor in your data structure
      m_point_cloud.emplace_back(point, normal, desc);
    }
  }
}

pcl::PointCloud<pcl::PointXYZ>::Ptr pointCloud::toPclCloud() {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
  cloud->reserve(m_point_cloud.size());
  for (size_t i = 0; i < m_point_cloud.size(); ++i) {
    pcl::PointXYZ pcl_pt;
    pcl_pt.x = m_point_cloud.at(i).s_point.x();
    pcl_pt.y = m_point_cloud.at(i).s_point.y();
    pcl_pt.z = m_point_cloud.at(i).s_point.z();

    cloud->push_back(pcl_pt);
  }

  return cloud;
}

pcl::PointCloud<pcl::PointXYZ>::Ptr
pointCloud::toPclCloud(const std::vector<adi::Point> &cloud) {
  pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_cloud(
      new pcl::PointCloud<pcl::PointXYZ>);
  int count = 0;
  std::cout << cloud.size() << std::endl;
  for (size_t i = 0; i < cloud.size(); ++i) {
    pcl::PointXYZ pcl_pt;
    if (static_cast<double>(std::isinf(cloud.at(i).s_point.x())) ||
        static_cast<double>(std::isnan(cloud.at(i).s_point.x())) ||
        static_cast<double>(std::isnan(cloud.at(i).s_point.y())) ||
        static_cast<double>(std::isinf(cloud.at(i).s_point.y())) ||
        static_cast<double>(std::isnan(cloud.at(i).s_point.z())) ||
        static_cast<double>(std::isinf(cloud.at(i).s_point.z()))) {
      count++;
      continue;
    }

    pcl_pt.x = cloud.at(i).s_point.x();
    pcl_pt.y = cloud.at(i).s_point.y();
    pcl_pt.z = cloud.at(i).s_point.z();

    pcl_cloud->push_back(pcl_pt);
  }

  std::cout << "Number of invalid points: " << count << std::endl;
  return pcl_cloud;
}

void pointCloud::loadPointCloud(
    const pcl::PointCloud<pcl::PointXYZ>::Ptr &point_cloud) {
  m_point_cloud.reserve(point_cloud->size());
  for (size_t i = 0; i < point_cloud->size(); ++i) {
    Point pt;
    pt.s_point.x() = point_cloud->at(i).x;
    pt.s_point.y() = point_cloud->at(i).y;
    pt.s_point.z() = point_cloud->at(i).z;
    m_point_cloud.push_back(pt);
  }
}

pcl::PointCloud<pcl::PointNormal>::Ptr
pointCloud::normalEstimation(pcl::PointCloud<pcl::PointXYZ>::Ptr point_cloud,
                             const double &search_radius) {
  if (point_cloud->points.size() == 0) {
    return nullptr;
  } else {
    pcl::PointCloud<pcl::Normal>::Ptr cloud_normals(
        new pcl::PointCloud<pcl::Normal>);
    pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> normalEstimation;
    normalEstimation.setNumberOfThreads(MAX_NUMBER_OF_THREADS);
    normalEstimation.setInputCloud(point_cloud);
    normalEstimation.setRadiusSearch(search_radius);

    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(
        new pcl::search::KdTree<pcl::PointXYZ>);
    normalEstimation.setSearchMethod(tree);

    pcl::PointXYZ centroid;
    pcl::computeCentroid<pcl::PointXYZ>(*point_cloud, centroid);
    normalEstimation.setViewPoint(centroid.x, centroid.y, centroid.z);

    pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
    normalEstimation.compute(*normals);

    for (auto &n : normals->points) {
      n.normal_x *= -1;
      n.normal_y *= -1;
      n.normal_z *= -1;
    }

    pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals(
        new pcl::PointCloud<pcl::PointNormal>);
    cloud_with_normals->reserve(point_cloud->size());
    for (size_t i = 0; i < point_cloud->size(); ++i) {
      // Skip if the normal contains NaN values
      if (std::isnan(normals->at(i).normal_x) ||
          std::isnan(normals->at(i).normal_y) ||
          std::isnan(normals->at(i).normal_z)) {
        continue;
      }

      pcl::PointNormal point_normal;
      point_normal.x = point_cloud->at(i).x;
      point_normal.y = point_cloud->at(i).y;
      point_normal.z = point_cloud->at(i).z;
      point_normal.normal_x = normals->at(i).normal_x;
      point_normal.normal_y = normals->at(i).normal_y;
      point_normal.normal_z = normals->at(i).normal_z;

      cloud_with_normals->push_back(point_normal);
    }
    return cloud_with_normals;
  }
}

void pointCloud::samplePointCloud(const uint32_t max_number_of_points) {
  if (max_number_of_points == 0) {
    return;
  }

  // Create a unique pointer for the subsampled point cloud
  std::unique_ptr<std::vector<adi::Point>> subsampled_point_cloud =
      std::make_unique<std::vector<adi::Point>>();
  subsampled_point_cloud->reserve(max_number_of_points);

  if (m_point_cloud.size() >= max_number_of_points) {
    // Initialize variables
    std::vector<bool> is_sampled(m_point_cloud.size(), false);
    std::vector<double> min_distances(m_point_cloud.size(),
                                      std::numeric_limits<double>::max());

    // Randomly pick the first point
    int random_index = std::rand() % m_point_cloud.size();
    subsampled_point_cloud->push_back(m_point_cloud[random_index]);
    is_sampled[random_index] = true;

    // Update distances for the first point
    for (size_t j = 0; j < m_point_cloud.size(); ++j) {
      if (j != random_index) {
        min_distances[j] =
            (m_point_cloud[j].s_point - m_point_cloud[random_index].s_point)
                .norm();
      }
    }

    // Iteratively sample the farthest points
    for (size_t i = 1; i < max_number_of_points; ++i) {
      // Find the farthest point using the min_distances array
      size_t farthest_index = 0;
      double max_distance = -std::numeric_limits<double>::max();
      for (size_t j = 0; j < m_point_cloud.size(); ++j) {
        if (!is_sampled[j] && min_distances[j] > max_distance) {
          max_distance = min_distances[j];
          farthest_index = j;
        }
      }

      // Add the farthest point to the sampled set
      subsampled_point_cloud->push_back(m_point_cloud[farthest_index]);
      is_sampled[farthest_index] = true;

      // Update the minimum distances for the remaining points
      const Eigen::Vector3d &new_sampled_point =
          m_point_cloud[farthest_index].s_point;
      for (size_t j = 0; j < m_point_cloud.size(); ++j) {
        if (!is_sampled[j]) {
          double distance =
              (m_point_cloud[j].s_point - new_sampled_point).norm();
          if (distance < min_distances[j]) {
            min_distances[j] = distance;
          }
        }
      }
    }
    this->setPointCloud(std::move(subsampled_point_cloud));
    return;
  }
}

void pointCloud::serializeCloud(const std::vector<adi::Point> &cloud,
                                const std::string &file_path) {
  pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_cloud = toPclCloud(cloud);

  pcl::io::savePLYFile(file_path, *pcl_cloud);
}
} // namespace adi
