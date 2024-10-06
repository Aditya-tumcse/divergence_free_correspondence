#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include <eigen3/Eigen/Dense>
#include <cmath>

#include "deformation_field.hpp"
#include "constants.hpp"

namespace adi{
    namespace numerics{
        const std::vector<adi::deformation_field::BasisIndices> GenerateBasisIndices(const uint32_t &max_number_of_velocity_basis);

        const double ObtainCoefficientOfVelocityField(const double eigen_value);

        const std::vector<std::map<double, adi::deformation_field::BasisIndices>> GenerateBaseIndexMap(const std::vector<adi::deformation_field::BasisIndices> &base_indices);

        Eigen::Vector3d RungeKutaIntegration(const Eigen::Vector3d &src_pt, const std::vector<adi::deformation_field::BasisIndices> &basis_indices,const Eigen::VectorXd &coeffs_ak,const uint32_t num_time_steps);
    
        Eigen::MatrixXd ComputeJacobian(const std::vector<adi::deformation_field::BasisIndices> &base_indices, const adi::Point &point);

        const double EucledianDistance(const Eigen::Vector3d &point_1, const Eigen::Vector3d &point_2);

        //TODO: Do not send in the copies of cloud_src and cloud_target
        Eigen::MatrixXd ComputeHessian(const std::vector<adi::deformation_field::BasisIndices> &base_indices, std::vector<adi::Point> cloud_src,std::vector<adi::Point> cloud_target, const Eigen::MatrixXd &soft_corr_matrix, const double r0);
    
        Eigen::VectorXd ComputeGradient(const std::vector<adi::deformation_field::BasisIndices> &base_indices, std::vector<adi::Point> cloud_src,std::vector<adi::Point> cloud_target, const Eigen::MatrixXd &soft_corr_matrix, const double r0);
    }
}
#endif