#ifndef UTILITIES_HPP
#define UTILTIES_HPP

#include "io.hpp"

namespace utilities{
    double eucledianDistance(const Eigen::Vector3d &point_1, const Eigen::Vector3d &point_2);

    Eigen::MatrixXd computeDistanceMatrix(const std::vector<Eigen::Vector3d> points);

    double computeL2Norm(const std::array<double, 352> &descriptor_1, const std::array<double, 352>  &descriptor_2);

    bool isValidDescriptor(const pcl::SHOT352& descriptor);

    std::unique_ptr<std::vector<std::pair<adi::Point, adi::Point>>> computeCorrespondences(std::vector<adi::Point> source_point_cloud,std::vector<adi::Point> target_point_cloud);

    
}
#endif