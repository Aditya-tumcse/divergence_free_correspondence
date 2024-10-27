#include "numerics.hpp"
#include "Fastor/Fastor.h"

#include <random>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

namespace adi{
 namespace numerics{
        const std::vector<adi::deformation_field::BasisIndices> GenerateBasisIndices(const uint32_t &max_number_of_velocity_basis)
        {
            std::vector<adi::deformation_field::BasisIndices> basis_indices;
            for(uint32_t i = 1;i <= max_number_of_velocity_basis;++i)
            {
                for(uint32_t j = 1;j <= max_number_of_velocity_basis;++j)
                {
                    for(uint32_t k = 1;k <= max_number_of_velocity_basis;++k)
                    {
                        adi::deformation_field::BasisIndices base_index;
                        base_index.index_1 = i;
                        base_index.index_2 = j;
                        base_index.index_3 = k;
                        double eigen_val = std::pow((std::pow(M_PI,2) * (std::pow(i,2) + std::pow(j,2) + std::pow(k,2))),-1.5);
                        base_index.eigen_value = eigen_val;
                        
                        basis_indices.push_back(base_index);
                    }
                }
            }

            return basis_indices;
        }

        const double ObtainCoefficientOfVelocityField(const double eigen_value)
        {
            const double mean = 0.0;
            const double stddev = eigen_value;

            // Create a random number generator and a normal distribution object
            std::random_device rd;  // Random device for seeding
            std::mt19937 gen(rd()); // Mersenne Twister generator seeded with rd()
            std::normal_distribution<> d(mean, stddev);

            // Generate a sample from the normal distribution
            double sample = d(gen);

            return sample;
        }

        // Eigen::Vector3d RungeKutaIntegration(const Eigen::Vector3d &src_pt, const std::vector<adi::deformation_field::BasisIndices> &basis_indices, const Eigen::VectorXd &coeffs_ak, const uint32_t num_time_steps)
        // {
        //     Eigen::Vector3d current_pt = src_pt;
        //     const double dt = 1.0 / num_time_steps;

        //     std::cout << "Size of each time step: " << dt << std::endl;

            
        //         Eigen::Vector3d local_current_pt = current_pt;
        //         for (uint32_t i = 0; i < num_time_steps; ++i)
        //         {
        //             adi::deformation_field::DeformationField df;

        //             Eigen::Vector3d k1 = df.computeVelocityField(local_current_pt, basis_indices, coeffs_ak);

        //             Eigen::Vector3d midpoint;
        //             midpoint.setZero();
        //             midpoint.x() = local_current_pt.x() + 0.5 * dt * k1.x();
        //             midpoint.y() = local_current_pt.y() + 0.5 * dt * k1.y();
        //             midpoint.z() = local_current_pt.z() + 0.5 * dt * k1.z();

        //             Eigen::Vector3d k2 = df.computeVelocityField(midpoint, basis_indices, coeffs_ak);

        //             current_pt.x() += dt * k2.x();
        //             current_pt.y() += dt * k2.y();
        //             current_pt.z() += dt * k2.z();
                    
        //         }

        //     return current_pt;
        // }

        Eigen::MatrixXd RungeKutta2Integration(const std::vector<adi::deformation_field::BasisIndices> &basis_indices, Eigen::MatrixXd src_point_cloud,
            const Eigen::VectorXd &coeffs_ak,
            const Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS, MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH> &vel_basis_functions,
            const uint32_t num_time_steps)
        {
            Eigen::MatrixXd updated_pts = src_point_cloud;

            const double dt = 1.0 / num_time_steps;
            
            std::cout << "Size of each time step: " << dt << std::endl;
            for(uint32_t t = 0;t < num_time_steps;++t)
            {
                adi::deformation_field::DeformationField df;
                Eigen::MatrixXd intermediate_vel_field = df.computeVelocityField(coeffs_ak, vel_basis_functions);

                Eigen::MatrixXd midpoint = Eigen::MatrixXd::Zero(NUMBER_OF_SAMPLE_POINTS,3);

                midpoint.col(0) = updated_pts.col(0) + 0.5 * dt * intermediate_vel_field.col(0);
                midpoint.col(1) = updated_pts.col(1) + 0.5 * dt * intermediate_vel_field.col(1);
                midpoint.col(2) = updated_pts.col(2) + 0.5 * dt * intermediate_vel_field.col(2);

                auto updated_vel_basis_functions = df.computeVelocityBasisFunctions(basis_indices, midpoint);

                Eigen::MatrixXd updated_vel_field = df.computeVelocityField(coeffs_ak, updated_vel_basis_functions);

                updated_pts.col(0) += dt * updated_vel_field.col(0);
                updated_pts.col(1) += dt * updated_vel_field.col(1);
                updated_pts.col(2) += dt * updated_vel_field.col(2);
            }
            
            return updated_pts; 
        }



        // Eigen::MatrixXd ComputeJacobian(const std::vector<adi::deformation_field::BasisIndices> &base_indices, const adi::Point &point)
        // {
        //     Eigen::MatrixXd jacobian(3,base_indices.size());
        //     adi::deformation_field::DeformationField df;

        //     for(uint32_t i = 0;i < base_indices.size();++i)
        //     {
        //         Eigen::Vector3d vel_basis_functions = df.computeVelocityBasisFunctions(base_indices[i].eigen_value, base_indices[i], point.s_point);
        //         jacobian(0, i) = vel_basis_functions.x();
        //         jacobian(1, i) = vel_basis_functions.y();
        //         jacobian(2, i) = vel_basis_functions.z();
        //     }

        //     return jacobian;
        // }

        const double EucledianDistance(const Eigen::Vector3d &point_1, const Eigen::Vector3d &point_2)
        {
            return (point_1 - point_2).norm();
        }

        // Eigen::MatrixXd ComputeHessian(const std::vector<adi::deformation_field::BasisIndices> &base_indices, std::vector<adi::Point> cloud_src,std::vector<adi::Point> cloud_target, const Eigen::MatrixXd &soft_corr_matrix, const double r0)
        // {
        //     Eigen::MatrixXd H = Eigen::MatrixXd::Zero(base_indices.size(), base_indices.size());

        //     for(uint32_t i = 0;i < cloud_src.size();++i)
        //     {
        //         Eigen::MatrixXd jacobian = ComputeJacobian(base_indices, cloud_src[i]);

        //         for(uint32_t j = 0;j < cloud_target.size();++j)
        //         {
        //             double eucledian_dist = EucledianDistance(cloud_src[i].s_point, cloud_target[j].s_point);
        //             double rho = (eucledian_dist <= r0) ? 0.5 * eucledian_dist * eucledian_dist : r0 * std::abs(eucledian_dist) - 0.5 * r0 * r0;

        //             H += soft_corr_matrix(i,j) * rho * (jacobian.transpose() * jacobian); 
        //         }
        //     }

        //     return H;
        // }

        // Eigen::VectorXd ComputeGradient(const std::vector<adi::deformation_field::BasisIndices> &base_indices, std::vector<adi::Point> cloud_src,std::vector<adi::Point> cloud_target, const Eigen::MatrixXd &soft_corr_matrix, const double r0)
        // {
        //     Eigen::VectorXd gradient = Eigen::VectorXd::Zero(base_indices.size());
        //     for(uint32_t i = 0;i < cloud_src.size();++i)
        //     {
        //         Eigen::MatrixXd jacobian = ComputeJacobian(base_indices, cloud_src[i]);

        //         for(uint32_t j = 0;j < cloud_target.size();++j)
        //         {
        //             double eucledian_dist = EucledianDistance(cloud_src[i].s_point, cloud_target[j].s_point);
        //             double rho = (eucledian_dist <= r0) ? 0.5 * eucledian_dist * eucledian_dist : r0 * std::abs(eucledian_dist) - 0.5 * r0 * r0;

        //             gradient += soft_corr_matrix(i,j) * rho * jacobian.transpose() * (cloud_src[i].s_point - cloud_target[j].s_point); 
        //         }
        //     }

        //     return gradient;
        // }


    }
}
