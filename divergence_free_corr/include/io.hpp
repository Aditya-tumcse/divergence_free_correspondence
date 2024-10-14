#ifndef IO_HPP
#define IO_HPP

#include<array>
#include <pcl/io/ply_io.h>
#include <pcl/point_cloud.h>
#include <pcl/common/io.h>
#include <pcl/point_types.h>
#include <pcl/features/shot_omp.h>
#include <pcl/features/normal_3d.h>
#include <Eigen/Dense>

namespace adi{

  struct Point{
      Eigen::Vector3d s_point; //Point Coordinates
      Eigen::Vector3d s_normal; // Normal info 
      std::array<double, 352> s_descriptor; // Feature descriptor

      Point() : s_point(0.0,0.0,0.0), s_normal(0.0,0.0,0.0), s_descriptor{} {}

      Point(Eigen::Vector3d point, Eigen::Vector3d normal, std::array<double, 352> descriptor) : s_point(point), s_normal(normal), s_descriptor(descriptor) {}

      Point(pcl::PointXYZ point) : s_normal(0.0,0.0,0.0), s_descriptor{}{
        s_point.x() = point.x;
        s_point.y() = point.y;
        s_point.z() = point.z;
      }

      pcl::SHOT352 convertToPCLDescriptor();
  };

  class pointCloud{
      public:
          pointCloud(const std::string &point_cloud_path);

          pointCloud(const pcl::PointCloud<pcl::PointXYZ>::Ptr &source_cloud);

          pcl::PointCloud<pcl::PointXYZ>::Ptr loadPlyFile(const std::string &cloud_file_path);

          void setPointCloud(const std::unique_ptr<std::vector<Point>> cloud);

          void computeShotFeatures(const double &search_radius);

          const std::vector<Point> getPointCloud() const{ return m_point_cloud;}

          void samplePointCloud(const uint32_t max_number_of_points);

          const uint32_t getSize() const {return m_point_cloud.size();}

          ~pointCloud() = default;

      private:
          std::vector<Point> m_point_cloud; 

          pcl::PointCloud<pcl::PointXYZ>::Ptr toPclCloud();
          
          void loadPointCloud(const pcl::PointCloud<pcl::PointXYZ>::Ptr &point_cloud);

          pcl::PointCloud<pcl::PointNormal>::Ptr normalEstimation(pcl::PointCloud<pcl::PointXYZ>::Ptr point_cloud, const double &search_radius);
          
  };
}
#endif // IO_HPP