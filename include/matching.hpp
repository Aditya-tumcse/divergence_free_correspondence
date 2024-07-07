#ifndef MATCHING_HPP
#define MATCHING_HPP

#include "io.hpp"

#include <eigen3/Eigen/Dense>

namespace adi{
    namespace matching{

        class Matching{
            public:
                Matching(const std::vector<std::pair<adi::Point, adi::Point>> &correspondences);

                Eigen::MatrixXd computeSoftCorrespondences(const std::vector<adi::Point> &source_point_cloud, const std::vector<adi::Point> &target_point_cloud);

            private:
                const std::vector<std::pair<adi::Point, adi::Point>> &m_correspondences;

                const double computeMeanEucledianDistance();
                const double computeMeanDescriptorDistance();
                const Eigen::MatrixXd computeMetricDistance(const std::vector<adi::Point> &source_point_cloud, const std::vector<adi::Point> &target_point_cloud);
        };
    }
}

#endif