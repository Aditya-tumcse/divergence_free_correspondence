#include "deformation_field.hpp"
#include "io.hpp"
#include "constants.hpp"

#include<cmath>

namespace adi{
    namespace deformation_field{
        deformationField::deformationField(const std::vector<adi::Point> &source_point_cloud, const std::vector<adi::Point> &target_point_cloud, const std::vector<std::pair<adi::Point, adi::Point>> &correspondences): m_source_point_cloud(source_point_cloud), m_target_point_cloud(target_point_cloud), m_correspondences(correspondences){}

        const double deformationField::computeMeanEucledianDistance()
        {
            double mean_eucledian_distance = 0;
            for(unsigned int i = 0;i < m_correspondences.size();++i)
            {
                mean_eucledian_distance += (m_correspondences[i].first.s_point - m_correspondences[i].second.s_point).norm();
            }
            return mean_eucledian_distance / m_correspondences.size();
        }

        const double deformationField::computeMeanDescriptorDistance()
        {
            double mean_descriptor_distance = 0;
            for(unsigned int i = 0;i < m_correspondences.size();++i)
            {
                mean_descriptor_distance += utilities::computeL2Norm(m_correspondences[i].first.s_descriptor, m_correspondences[i].second.s_descriptor);
            }
            return mean_descriptor_distance / m_correspondences.size();
        }

        const Eigen::MatrixXd deformationField::computeMetricDistance()
        {
            Eigen::MatrixXd metric_distance = Eigen::MatrixXd::Zero(m_source_point_cloud.size(), m_source_point_cloud.size());
            const double mean_eucledian_distance = this->computeMeanEucledianDistance();
            const double mean_descriptor_distance = this->computeMeanDescriptorDistance();
            const double factor = mean_eucledian_distance / mean_descriptor_distance;

            for(unsigned int i = 0;i < m_source_point_cloud.size();++i)
            {
                for(unsigned int j = 0;j < m_target_point_cloud.size();++j)
                {
                    const double local_eucledian_distance = (m_source_point_cloud[i].s_point - m_target_point_cloud[j].s_point).norm();
                    const double local_descriptor_distance = utilities::computeL2Norm(m_source_point_cloud[i].s_descriptor, m_target_point_cloud[j].s_descriptor);
                    metric_distance(i,j) = local_eucledian_distance + factor * local_descriptor_distance;
                }
            }

            return metric_distance;
        }
        
        Eigen::MatrixXd deformationField::computeSoftCorrespondences()
        {
            Eigen::MatrixXd correspondences = Eigen::MatrixXd::Zero(m_source_point_cloud.size(), m_source_point_cloud.size());
            const Eigen::MatrixXd metric_distance = computeMetricDistance();

            for(unsigned int i = 0;i < m_source_point_cloud.size();++i)
            {
                for(unsigned int j = 0;j < m_target_point_cloud.size();++j)
                {
                    const double distance = metric_distance(i,j);

                    const double numerator = std::exp((-1/(2 * SIGMA * SIGMA) * distance * distance));
                    const double denominator = std::exp(metric_distance.row(i).array().sum() * (-1 / (2 * SIGMA * SIGMA))) + std::pow((2 * M_PI * SIGMA * SIGMA), 1.5);
                    correspondences(i,j) = numerator / denominator;
                }
            }
            return correspondences;
        }
    }
}