#include "deformation_field.hpp"
#include "io.hpp"
#include "constants.hpp"

#include<cmath>
#include<random>

namespace adi{
    namespace deformation_field{
        const double DeformationField::phi(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z )
        {
            return 0.5 * sin(M_PI * index_1 * x) * sin(M_PI * index_2 * y) * sin(M_PI * index_3 * z);
        }

        const double DeformationField::dphidx(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z )
        {
            return 0.5 * M_PI * index_1 * cos(M_PI * index_1 * x) * sin(M_PI * index_2 * y) * sin(M_PI * index_3 * z);
        }

        const double DeformationField::dphidy(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z )
        {
            return 0.5 * M_PI * index_2 * cos(M_PI * index_1 * x) * sin(M_PI * index_2 * y) * sin(M_PI * index_3 * z);
        }

        const double DeformationField::dphidz(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z )
        {
            return 0.5 * M_PI * index_3 * cos(M_PI * index_1 * x) * sin(M_PI * index_2 * y) * sin(M_PI * index_3 * z);
        }

        const std::vector<BasisIndices> DeformationField::GenerateBasisIndices(const uint32_t &max_number_of_velocity_basis)
        {
            std::vector<BasisIndices> basis_indices;
            for(uint32_t i = 1;i <= max_number_of_velocity_basis;++i)
            {
                for(uint32_t j = 1;j <= max_number_of_velocity_basis;++j)
                {
                    for(uint32_t k = 1;k <= max_number_of_velocity_basis;++k)
                    {
                        BasisIndices base_index;
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

        const double DeformationField::ObtainCoefficientOfVelocityField(const double eigen_value)
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

        std::vector<Eigen::Vector3d> DeformationField::computeVelocityField()
        {
            std::vector<Eigen::Vector3d> velocity_field;
            const double grid_spacing = DOMAIN_LENGTH / (m_number_of_grid_points - 1);
            std::vector<BasisIndices> base_indices = GenerateBasisIndices(MAX_NUMBER_OF_VELOCITY_BASIS);
            for(uint32_t i = 0;i < m_number_of_grid_points;++i)
            {
                double x = i * grid_spacing;
                for(uint32_t j = 0;j < m_number_of_grid_points;++j)
                {
                    double y = j * grid_spacing;
                    for(uint32_t k = 0;k<m_number_of_grid_points;++k)
                    {
                        double z = k * grid_spacing;

                        // Initialize velocity components
                        double vx_sum = 0.0, vy_sum = 0.0, vz_sum = 0.0;
                        Eigen::Vector3d velocity_field_per_point;
                        for(uint32_t b = 0;b < base_indices.size();++b)
                        {
                            auto base_index = base_indices[b];

                            // Compute partial derivatives of the eigenfunction phi
                            double dphi_dx = dphidx(base_index.index_1, base_index.index_2, base_index.index_3, x, y, z);
                            double dphi_dy = dphidy(base_index.index_1, base_index.index_2, base_index.index_3, x, y, z);
                            double dphi_dz = dphidz(base_index.index_1, base_index.index_2, base_index.index_3, x, y, z);

                            // Basis function velocities
                            double v1x = 0, v1y = dphi_dz, v1z = -dphi_dy;
                            double v2x = -dphi_dz, v2y = 0, v2z = dphi_dx;
                            double v3x = dphi_dy, v3y = -dphi_dx, v3z = 0;

                            const double coeff_ak = ObtainCoefficientOfVelocityField(base_index.eigen_value);

                            vx_sum += coeff_ak * (v1x + v2x + v3x);
                            vy_sum += coeff_ak * (v1y + v2y + v3y);
                            vz_sum += coeff_ak * (v1z + v2z + v3z);
                        }
                        velocity_field_per_point.x() = vx_sum;
                        velocity_field_per_point.y() = vy_sum;
                        velocity_field_per_point.z() = vz_sum;

                        velocity_field.emplace_back(velocity_field_per_point);
                    }
                }
            }
            return velocity_field;
        }

    }
}