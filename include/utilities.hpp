#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/registration/correspondence_estimation.h>

namespace utilities{

    double eucledianDistance(const Eigen::Vector3d &point_1, const Eigen::Vector3d &point_2)
    {
        return (std::sqrt(std::pow((point_1[0] - point_2[0]), 2) + std::pow((point_1[1] - point_2[1]), 2) + std::pow((point_1[2] - point_2[2]), 2)));
    }

    Eigen::MatrixXd computeDistanceMatrix(const std::vector<Eigen::Vector3d> points){
        Eigen::MatrixXd distance_matrix = Eigen::MatrixXd::Zero(points.size(), points.size());
        for(int i = 0;i < points.size();++i){
            for(int j = 0;j < points.size();++j){
                double eucledian_distance = eucledianDistance(points[i], points[j]);
                distance_matrix(i,j) = eucledian_distance;
            }
        }
        return distance_matrix;
    }

    std::vector<std::pair<adi::Point, adi::Point>> computeCorrespondences(std::vector<adi::Point> &source_point_cloud,std::vector<adi::Point> &target_point_cloud)
    {

        assert(source_point_cloud.size() == target_point_cloud.size());
        std::vector<std::pair<adi::Point, adi::Point>> correspondences;
        pcl::PointCloud<pcl::SHOT352>::Ptr target_descriptors(new pcl::PointCloud<pcl::SHOT352>);
        for(adi::Point &point : target_point_cloud)
        {
            target_descriptors->push_back(point.convertToPCLDescriptor());
        }
        pcl::KdTreeFLANN<pcl::SHOT352> kdtree;
        kdtree.setInputCloud(target_descriptors);

        // Find correspondence
        for(auto &source_point : source_point_cloud)
        {
            pcl::SHOT352 source_descriptor = source_point.convertToPCLDescriptor();
            std::vector<int> nn_indices(1);
            std::vector<float> nn_distances(1);

            if(kdtree.nearestKSearch(source_descriptor, 1, nn_indices, nn_distances) > 0)
            {
                int target_index = nn_indices[0];
                correspondences.emplace_back(std::make_pair(source_point, target_point_cloud[target_index]));
            }
        }
        return correspondences;
    }


}

#endif