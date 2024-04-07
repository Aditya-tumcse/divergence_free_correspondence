#ifndef IO_HPP
#define IO_HPP

#include<array>
#include <pcl/io/ply_io.h>
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
};


class pointCloud{
    public:
        pointCloud(std::string point_cloud_path) : m_point_cloud_path(point_cloud_path){}

        const std::string getPointCloudPath() const {return m_point_cloud_path;}

        pcl::PointCloud<pcl::PointXYZ>::Ptr readPointCloud(const std::string &point_cloud_path);

        pcl::PointCloud<pcl::PointNormal>::Ptr normalEstimation(pcl::PointCloud<pcl::PointXYZ>::Ptr point_cloud, const float &search_radius);

        void computeShotFeatures(pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals);

        const std::vector<Point> getPointCloud() const{ return m_point_cloud;}

        std::vector<Point> samplePointCloud(std::vector<adi::Point> point_cloud, const int max_number_of_points);

    private:
        std::string m_point_cloud_path;
        std::vector<Point> m_point_cloud; 

        const std::vector<Eigen::Vector3d> extractPoints(std::vector<adi::Point> points);
};
}
#endif // IO_HPP