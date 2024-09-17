#ifndef IO_HPP
#define IO_HPP

#include<array>
#include <pcl/io/ply_io.h>
#include <pcl/common/io.h>
#include <pcl/point_types.h>
#include <pcl/features/shot_omp.h>
#include <pcl/features/normal_3d.h>
#include <Eigen/Dense>
#include <queue>

namespace adi{

struct Point{
    Eigen::Vector3d s_point; //Point Coordinates
    Eigen::Vector3d s_normal; // Normal info 
    std::array<double, 352> s_descriptor; // Feature descriptor

    Point() : s_point(0.0,0.0,0.0), s_normal(0.0,0.0,0.0), s_descriptor{} {}

    Point(Eigen::Vector3d point, Eigen::Vector3d normal, std::array<double, 352> descriptor) : s_point(point), s_normal(normal), s_descriptor(descriptor) {}

    pcl::SHOT352 convertToPCLDescriptor();
};

// Comparator for priority queue (max heap)
struct FarthestPointComparator {
  bool operator()(const std::pair<double, size_t> &a, const std::pair<double, size_t> &b) {
    return a.first < b.first; // max heap based on distance
  }
};

class pointCloud{
    public:
        pointCloud(const std::string &point_cloud_path, const double &search_radius);

        const std::vector<Point> getPointCloud() const{ return m_point_cloud;}

        uint32_t getColumnIdOfTheFarthestSample(const Eigen::RowVectorXd &row);
        std::vector<Point> samplePointCloud(const uint32_t max_number_of_points);

        const uint32_t getSize() const {return m_point_cloud.size();}

    private:
        std::vector<Point> m_point_cloud; 
        
        pcl::PointCloud<pcl::PointXYZ>::Ptr readPointCloud(const std::string &point_cloud_path);
        pcl::PointCloud<pcl::PointNormal>::Ptr normalEstimation(pcl::PointCloud<pcl::PointXYZ>::Ptr point_cloud, const double &search_radius);
        void computeShotFeatures(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_with_normals, const double &search_radius);
        
        const std::vector<Eigen::Vector3d> extractPoints(std::vector<adi::Point> points);
};
}

namespace utilities{
    double eucledianDistance(const Eigen::Vector3d &point_1, const Eigen::Vector3d &point_2);

    Eigen::MatrixXd computeDistanceMatrix(const std::vector<Eigen::Vector3d> points);

    double computeL2Norm(const std::array<double, 352> &descriptor_1, const std::array<double, 352>  &descriptor_2);

    bool isValidDescriptor(const pcl::SHOT352& descriptor);

    std::vector<std::pair<adi::Point, adi::Point>> computeCorrespondences(std::vector<adi::Point> source_point_cloud,std::vector<adi::Point> target_point_cloud);

    
}
#endif // IO_HPP