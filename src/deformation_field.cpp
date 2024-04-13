#include "deformation_field.hpp"
#include "io.hpp"
#include "constants.hpp"

#include<cmath>

namespace adi{
    namespace deformation_field{

        const std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> getBasisFunctions(unsigned int number_of_basis_functions)
        {
             // Define a lambda function to calculate the sum of squares of tuple elements
            auto sum_of_squares = [](const std::tuple<unsigned int, unsigned int, unsigned int>& t) {
                unsigned int sum = 0;
                sum += std::get<0>(t) * std::get<0>(t);
                sum += std::get<1>(t) * std::get<1>(t);
                sum += std::get<2>(t) * std::get<2>(t);
                return sum;
            };

            std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> basis_functions;
            const unsigned int number_of_basis_functions_per_dimension = number_of_basis_functions / 3;
            for(unsigned int i = 0;i < number_of_basis_functions_per_dimension / 3;++i)
            {
                for(unsigned int j = 0;j < number_of_basis_functions_per_dimension;++j)
                {
                    for(unsigned int k = 0;k < number_of_basis_functions_per_dimension;++k)
                    {
                        basis_functions.emplace_back(std::make_tuple(i,j,k));
                    }
                }
            }

            std::sort(basis_functions.begin(), basis_functions.end(), [&](const auto& a, const auto& b) {
            return sum_of_squares(a) > sum_of_squares(b);});

            // Take the top 1000 elements
            std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> top_1000(basis_functions.begin(), basis_functions.begin() + 1000);

            return top_1000;
        }

        deformationField::deformationField(const std::vector<adi::Point> &source_point_cloud, const std::vector<adi::Point> &target_point_cloud, const std::vector<std::pair<adi::Point, adi::Point>> &correspondences): m_source_point_cloud(source_point_cloud), m_target_point_cloud(target_point_cloud), m_correspondences(correspondences){}
        
        void deformationField::initDeformationField(const std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> &indices_for_basis_function_computation, const std::vector<adi::Point> &point_cloud)
        {
            Eigen::Matrix3d gradient_of_basis_functions = this->computeBasisFunctions(point_cloud, indices_for_basis_function_computation);
        }
        
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

        const Eigen::Matrix3d deformationField::computeBasisFunctions(const std::vector<adi::Point> &point_cloud, const std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> &indices_for_basis_function_computation)
        {
            Eigen::Matrix3d basis_functions = Eigen::Matrix3d::Zero(point_cloud.size(), 3);
            for(unsigned int i = 0;i < point_cloud.size();++i)
            {
                basis_functions(i,0) = static_cast<double>(0.5 * std::sin(point_cloud[i].s_point[0] * M_PI * std::get<0>(indices_for_basis_function_computation[i])));
                basis_functions(i,1) = static_cast<double>(0.5 * std::sin(point_cloud[i].s_point[1] * M_PI * std::get<1>(indices_for_basis_function_computation[i])));
                basis_functions(i,2) = static_cast<double>(0.5 * std::sin(point_cloud[i].s_point[2] * M_PI * std::get<2>(indices_for_basis_function_computation[i])));
            }
            return basis_functions;
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