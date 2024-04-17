#include "deformation_field.hpp"
#include "io.hpp"
#include "constants.hpp"

#include<cmath>

namespace adi{
    namespace deformation_field{

        const std::vector<Indices> getBasisFunctions(unsigned int number_of_basis_functions)
        {
             // Define a lambda function to calculate the sum of squares of tuple elements
            auto sum_of_squares = [](const Indices& t) {
                unsigned int sum = 0;
                sum += t.s_x * t.s_x;
                sum += t.s_y * t.s_y;
                sum += t.s_z * t.s_z;
                return sum;
            };

            std::vector<Indices> basis_functions;
            const unsigned int number_of_basis_functions_per_dimension = number_of_basis_functions / 3;
            for(unsigned int i = 0;i < number_of_basis_functions_per_dimension / 3;++i)
            {
                for(unsigned int j = 0;j < number_of_basis_functions_per_dimension;++j)
                {
                    for(unsigned int k = 0;k < number_of_basis_functions_per_dimension;++k)
                    {
                        basis_functions.emplace_back(Indices(i,j,k));
                    }
                }
            }

            std::sort(basis_functions.begin(), basis_functions.end(), [&](const auto& a, const auto& b) {
            return sum_of_squares(a) > sum_of_squares(b);});

            // Take the top 1000 elements
            std::vector<adi::deformation_field::Indices> top_1000(basis_functions.begin(), basis_functions.begin() + 1000);

            return top_1000;
        }

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

        const std::vector<double> deformationField::computGradientOfBasisFunctions(const std::vector<adi::Point> &point_cloud, const std::vector<Indices> &indices_for_basis_function_computation, const unsigned int &dimension)
        {
            std::vector<double> gradient_of_basis_functions;
            gradient_of_basis_functions.reserve(point_cloud.size());

            const unsigned int k = dimension;
            const unsigned int l = (k + 1) % 3;
            const unsigned int m = (k + 2) % 3;
            for(unsigned int i = 0;i < point_cloud.size();++i)
            {
                for(unsigned int j = 0;j < indices_for_basis_function_computation.size();++j)
                {
                    double result;
                    if(dimension %3 == 0)
                        result = std::pow(0.5,3) * M_PI * indices_for_basis_function_computation[j].s_x * std::cos(M_PI * point_cloud[i].s_point[k] * indices_for_basis_function_computation[j].s_x) * std::sin(M_PI * point_cloud[i].s_point[l] * indices_for_basis_function_computation[j].s_y) * std::sin(M_PI * point_cloud[i].s_point[m] * indices_for_basis_function_computation[j].s_z);
                    else if (dimension %3 == 1)
                        result = std::pow(0.5,3) * M_PI * indices_for_basis_function_computation[j].s_y * std::cos(M_PI * point_cloud[i].s_point[k] * indices_for_basis_function_computation[j].s_y) * std::sin(M_PI * point_cloud[i].s_point[l] * indices_for_basis_function_computation[j].s_x) * std::sin(M_PI * point_cloud[i].s_point[m] * indices_for_basis_function_computation[j].s_z);
                    else
                        result = std::pow(0.5,3) * M_PI * indices_for_basis_function_computation[j].s_z * std::cos(M_PI * point_cloud[i].s_point[k] * indices_for_basis_function_computation[j].s_z) * std::sin(M_PI * point_cloud[i].s_point[l] * indices_for_basis_function_computation[j].s_x) * std::sin(M_PI * point_cloud[i].s_point[m] * indices_for_basis_function_computation[j].s_y);
                    gradient_of_basis_functions.emplace_back(result);
                }
            }
            return gradient_of_basis_functions;
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