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

        const double EucledianDistance(const Eigen::Vector3d &point_1, const Eigen::Vector3d &point_2)
        {
            return (point_1 - point_2).norm();
        }

        struct GaussianPriorCostFunctor {
            GaussianPriorCostFunctor(const Eigen::MatrixXd &L_inv)
                :  s_L_inv(L_inv) {}

            template <typename T>
            bool operator()(const T* coeffs_ak, T* residual) const {
                // Compute the Gaussian prior residual term with (1/2) * (a^T * L_inv * a)
                Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> coeffs(coeffs_ak, s_L_inv.rows());
                Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> s_L_inv_T = s_L_inv.template cast<T>();

                T quadratic_term = coeffs.transpose() * s_L_inv_T * coeffs; // Ensure L_inv is cast to type T
                residual[0] = T(0.5) * quadratic_term; 

                return true;
            }

            const Eigen::MatrixXd &s_L_inv; // Reference to the L_inv matrix
        };

        void Optimize(const Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS, MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH> &vel_basis_functions, Eigen::VectorXd *coeffs_ak, const std::vector<adi::Point> &target_cloud, const Eigen::MatrixXd &soft_corr_matrix,const Eigen::MatrixXd& L_inv)
        {
            // Prepare the ceres problem
            ceres::Problem problem;

            // for(size_t i = 0;i < target_cloud.size();++i)
            // {
            //     auto weight = soft_corr_matrix.col(i);

            //     // Create a residual block for this target point
            //     problem.AddResidualBlock(
            //         new ceres::AutoDiffCostFunction<MatchingCostFunctor, 3, Eigen::Dynamic>(new MatchingCostFunctor(target_cloud[i], vel_basis_functions, coeffs_ak, weight)),
            //         nullptr,
            //         coeffs_ak->data() 
            //     );
            // }

            // Add Gaussian prior residual block
            problem.AddResidualBlock(
                new ceres::AutoDiffCostFunction<GaussianPriorCostFunctor, 1, MAX_NUMBER_OF_VELOCITY_BASIS>(new GaussianPriorCostFunctor(L_inv)),
                nullptr, 
                coeffs_ak->data() // Optimizing coeffs_ak
            );

            // Set solver options
            ceres::Solver::Options options;
            options.max_num_iterations = 10;
            options.linear_solver_type = ceres::DENSE_QR;

            // Solve the problem
            ceres::Solver::Summary summary;
            ceres::Solve(options, &problem, &summary);
            
            // Print solver summary
            std::cout << summary.FullReport() << std::endl;
        }

    }
}
