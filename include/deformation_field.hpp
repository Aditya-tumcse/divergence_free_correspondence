#ifndef DEFORMATION_FIELD_HPP
#define DEFORMATION_FIELD_HPP

#include "io.hpp"

#include <eigen3/Eigen/Dense>
#include <cassert>

namespace adi{
    namespace deformation_field{
        class deformationField{
            public:
                deformationField(const std::vector<adi::Point> &source_point_cloud, const std::vector<adi::Point> &target_point_cloud, const std::vector<std::pair<adi::Point, adi::Point>> &correspondences);

                Eigen::MatrixXd computeSoftCorrespondences();

            private:
                const std::vector<adi::Point> &m_source_point_cloud;
                const std::vector<adi::Point> &m_target_point_cloud;
                const std::vector<std::pair<adi::Point, adi::Point>> &m_correspondences;

                const double computeMeanEucledianDistance();
                const double computeMeanDescriptorDistance();
                const Eigen::MatrixXd computeMetricDistance();
        };
    }
}
#endif
