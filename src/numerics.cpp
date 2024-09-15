#include <random>

#include "numerics.hpp"

namespace adi{
 namespace numerics{
        const std::vector<adi::deformation_field::BasisIndices> GenerateBasisIndices(const uint32_t &max_number_of_velocity_basis)
        {
            std::vector<adi::deformation_field::BasisIndices> basis_indices;
            for(uint32_t i = 1;i <= max_number_of_velocity_basis / 3;++i)
            {
                for(uint32_t j = 1;j <= max_number_of_velocity_basis / 3;++j)
                {
                    for(uint32_t k = 1;k <= max_number_of_velocity_basis / 3;++k)
                    {
                        adi::deformation_field::BasisIndices base_index;
                        base_index.index_1 = i;
                        base_index.index_2 = j;
                        base_index.index_3 = k;
                        double eigen_val = -std::pow(M_PI,2) * (std::pow(i,2) + std::pow(j,2) + std::pow(k,2));
                        base_index.eigen_value = eigen_val;
                        
                        basis_indices.emplace_back(base_index);
                    }
                }
            }

            return basis_indices;
        }

        //TODO: Figure why this is being done
        const std::vector<std::map<double, adi::deformation_field::BasisIndices>> GenerateBaseIndexMap(const std::vector<adi::deformation_field::BasisIndices> &base_indices)
        {
            std::vector<std::map<double, adi::deformation_field::BasisIndices>> base_index_map;
            for(uint32_t i = 0;i < base_indices.size();++i)
            {
                std::map<double, adi::deformation_field::BasisIndices> base_id_map;
                base_id_map.insert(std::make_pair(0,base_indices.at(i)));
            }
            return base_index_map;
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

        Eigen::Vector3d RungeKutaIntegration(const Eigen::Vector3d &src_pt, const std::vector<adi::deformation_field::BasisIndices> &basis_indices,const Eigen::VectorXd &coeffs_ak,const double dt)
        {
            Eigen::Vector3d updated_pt;
            updated_pt.setZero();

            adi::deformation_field::DeformationField df;

            Eigen::Vector3d k1 = df.computeVelocityField(src_pt, basis_indices, coeffs_ak);

            Eigen::Vector3d midpoint;
            midpoint.setZero();
            midpoint.x() = src_pt.x() + 0.5 * dt * k1.x();
            midpoint.y() = src_pt.y() + 0.5 * dt * k1.y();
            midpoint.z() = src_pt.z() + 0.5 * dt * k1.z();

            Eigen::Vector3d k2 = df.computeVelocityField(midpoint, basis_indices, coeffs_ak);
            
            updated_pt.x() = src_pt.x() + dt * k2.x();
            updated_pt.y() = src_pt.y() + dt * k2.y();
            updated_pt.z() = src_pt.z() + dt * k2.z();

            return updated_pt;
        }

        Eigen::MatrixXd ComputeJacobian(const std::vector<adi::deformation_field::BasisIndices> &base_indices, const adi::Point &point)
        {
            Eigen::MatrixXd jacobian(3,base_indices.size());
            adi::deformation_field::DeformationField df;

            for(uint32_t i = 0;i < base_indices.size();++i)
            {
                std::tuple<double, double, double> vel_basis_functions = df.computeVelocityBasisFunctions(base_indices[i].eigen_value, base_indices[i], point.s_point);
                jacobian(0, i) = std::get<0>(vel_basis_functions);
                jacobian(1, i) = std::get<1>(vel_basis_functions);
                jacobian(2, i) = std::get<2>(vel_basis_functions);
            }

            return jacobian;
        }

        const double EucledianDistance(const Eigen::Vector3d &point_1, const Eigen::Vector3d &point_2)
        {
            return (point_1 - point_2).norm();
        }

        Eigen::MatrixXd ComputeHessian(const std::vector<adi::deformation_field::BasisIndices> &base_indices, std::vector<adi::Point> cloud_src,std::vector<adi::Point> cloud_target, const Eigen::MatrixXd &soft_corr_matrix, const double r0)
        {
            Eigen::MatrixXd H = Eigen::MatrixXd::Zero(base_indices.size(), base_indices.size());

            for(uint32_t i = 0;i < cloud_src.size();++i)
            {
                Eigen::MatrixXd jacobian = ComputeJacobian(base_indices, cloud_src[i]);

                for(uint32_t j = 0;j < cloud_target.size();++j)
                {
                    double eucledian_dist = EucledianDistance(cloud_src[i].s_point, cloud_target[j].s_point);
                    double rho = (eucledian_dist <= r0) ? 0.5 * eucledian_dist * eucledian_dist : r0 * std::abs(eucledian_dist) - 0.5 * r0 * r0;

                    H += soft_corr_matrix(i,j) * rho * (jacobian.transpose() * jacobian); 
                }
            }

            return H;
        }

        Eigen::VectorXd ComputeGradient(const std::vector<adi::deformation_field::BasisIndices> &base_indices, std::vector<adi::Point> cloud_src,std::vector<adi::Point> cloud_target, const Eigen::MatrixXd &soft_corr_matrix, const double r0)
        {
            Eigen::VectorXd gradient = Eigen::VectorXd::Zero(base_indices.size());
            for(uint32_t i = 0;i < cloud_src.size();++i)
            {
                Eigen::MatrixXd jacobian = ComputeJacobian(base_indices, cloud_src[i]);

                for(uint32_t j = 0;j < cloud_target.size();++j)
                {
                    double eucledian_dist = EucledianDistance(cloud_src[i].s_point, cloud_target[j].s_point);
                    double rho = (eucledian_dist <= r0) ? 0.5 * eucledian_dist * eucledian_dist : r0 * std::abs(eucledian_dist) - 0.5 * r0 * r0;

                    gradient += soft_corr_matrix(i,j) * rho * jacobian.transpose() * (cloud_src[i].s_point - cloud_target[j].s_point); 
                }
            }

            return gradient;
        }


    }
}
