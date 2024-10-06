#include "matching.hpp"
#include "constants.hpp"

namespace adi{
    namespace matching{

        Matching::Matching(const std::vector<std::pair<adi::Point, adi::Point>> &correspondences) : m_correspondences(correspondences){}
        
        const double Matching::computeMeanEucledianDistance()
        {
            double mean_eucledian_distance = 0;
            for(uint32_t i = 0;i < m_correspondences.size();++i)
            {
                mean_eucledian_distance += (m_correspondences[i].first.s_point - m_correspondences[i].second.s_point).norm();
            }
            return mean_eucledian_distance / m_correspondences.size();
        }

        const double Matching::computeMeanDescriptorDistance()
        {
            double mean_descriptor_distance = 0;
            for(uint32_t i = 0;i < m_correspondences.size();++i)
            {
                mean_descriptor_distance += utilities::computeL2Norm(m_correspondences[i].first.s_descriptor, m_correspondences[i].second.s_descriptor);
            }
            return mean_descriptor_distance / m_correspondences.size();
        }

        const Eigen::MatrixXd Matching::computeMetricDistance(const std::vector<adi::Point> &source_point_cloud, const std::vector<adi::Point> &target_point_cloud)
        {
            Eigen::MatrixXd metric_distance = Eigen::MatrixXd::Zero(source_point_cloud.size(), source_point_cloud.size());
            const double mean_eucledian_distance = this->computeMeanEucledianDistance();
            const double mean_descriptor_distance = this->computeMeanDescriptorDistance();
            const double factor = mean_eucledian_distance / mean_descriptor_distance;

            for(uint32_t i = 0;i < source_point_cloud.size();++i)
            {
                for(uint32_t j = 0;j < target_point_cloud.size();++j)
                {
                    const double local_eucledian_distance = (source_point_cloud[i].s_point - target_point_cloud[j].s_point).norm();
                    const double local_descriptor_distance = utilities::computeL2Norm(source_point_cloud[i].s_descriptor, target_point_cloud[j].s_descriptor);
                    metric_distance(i,j) = local_eucledian_distance + factor * local_descriptor_distance;
                }
            }

            return metric_distance;
        }
        
        Eigen::MatrixXd Matching::computeSoftCorrespondences(const std::vector<adi::Point> &source_point_cloud, const std::vector<adi::Point> &target_point_cloud)
        {
            Eigen::MatrixXd correspondences = Eigen::MatrixXd::Zero(source_point_cloud.size(), source_point_cloud.size());
            const Eigen::MatrixXd metric_distance = computeMetricDistance(source_point_cloud, target_point_cloud);

            for(uint32_t i = 0;i < source_point_cloud.size();++i)
            {
                for(uint32_t j = 0;j < target_point_cloud.size();++j)
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