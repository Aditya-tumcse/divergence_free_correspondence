#include "deformation_field.hpp"
#include "io.hpp"
#include "constants.hpp"
#include "numerics.hpp"

#include<cmath>


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

        // std::vector<Eigen::Vector3d> DeformationField::computeVelocityField(const Eigen::Vector3d &point)
        // {
        //     std::vector<Eigen::Vector3d> velocity_field;
        //     const double grid_spacing = DOMAIN_LENGTH / (m_number_of_grid_points - 1);
        //     std::vector<BasisIndices> base_indices = GenerateBasisIndices(MAX_NUMBER_OF_VELOCITY_BASIS);
        //     for(uint32_t i = 0;i < m_number_of_grid_points;++i)
        //     {
        //         double x = i * grid_spacing;
        //         for(uint32_t j = 0;j < m_number_of_grid_points;++j)
        //         {
        //             double y = j * grid_spacing;
        //             for(uint32_t k = 0;k<m_number_of_grid_points;++k)
        //             {
        //                 double z = k * grid_spacing;

        //                 // Initialize velocity components
        //                 double vx_sum = 0.0, vy_sum = 0.0, vz_sum = 0.0;
        //                 Eigen::Vector3d velocity_field_per_point;
        //                 for(uint32_t b = 0;b < base_indices.size();++b)
        //                 {
        //                     auto base_index = base_indices[b];

        //                     // Compute partial derivatives of the eigenfunction phi
        //                     double dphi_dx = dphidx(base_index.index_1, base_index.index_2, base_index.index_3, x, y, z);
        //                     double dphi_dy = dphidy(base_index.index_1, base_index.index_2, base_index.index_3, x, y, z);
        //                     double dphi_dz = dphidz(base_index.index_1, base_index.index_2, base_index.index_3, x, y, z);

        //                     // Basis function velocities
        //                     double v1x = 0, v1y = dphi_dz, v1z = -dphi_dy;
        //                     double v2x = -dphi_dz, v2y = 0, v2z = dphi_dx;
        //                     double v3x = dphi_dy, v3y = -dphi_dx, v3z = 0;

        //                     const double coeff_ak = ObtainCoefficientOfVelocityField(base_index.eigen_value);

        //                     vx_sum += coeff_ak * (v1x + v2x + v3x);
        //                     vy_sum += coeff_ak * (v1y + v2y + v3y);
        //                     vz_sum += coeff_ak * (v1z + v2z + v3z);
        //                 }
        //                 velocity_field_per_point.x() = vx_sum;
        //                 velocity_field_per_point.y() = vy_sum;
        //                 velocity_field_per_point.z() = vz_sum;

        //                 velocity_field.emplace_back(velocity_field_per_point);
        //             }
        //         }
        //     }
        //     return velocity_field;
        // }

        Eigen::Vector3d DeformationField::computeVelocityField(const Eigen::Vector3d &point, const std::vector<BasisIndices> &base_indices)
        {
            Eigen::Vector3d velocity_field;

            // Initialize velocity components
            double vx_sum = 0.0, vy_sum = 0.0, vz_sum = 0.0;
    
            for(uint32_t b = 0;b < base_indices.size();++b)
            {
                auto base_index = base_indices[b];

                // Compute partial derivatives of the eigenfunction phi
                double dphi_dx = dphidx(base_index.index_1, base_index.index_2, base_index.index_3, point.x(), point.y(), point.z());
                double dphi_dy = dphidy(base_index.index_1, base_index.index_2, base_index.index_3, point.x(), point.y(), point.z());
                double dphi_dz = dphidz(base_index.index_1, base_index.index_2, base_index.index_3, point.x(), point.y(), point.z());

                // Basis function velocities
                double v1x = 0, v1y = dphi_dz, v1z = -dphi_dy;
                double v2x = -dphi_dz, v2y = 0, v2z = dphi_dx;
                double v3x = dphi_dy, v3y = -dphi_dx, v3z = 0;

                const double coeff_ak = numerics::ObtainCoefficientOfVelocityField(base_index.eigen_value);

                vx_sum += coeff_ak * (v1x + v2x + v3x);
                vy_sum += coeff_ak * (v1y + v2y + v3y);
                vz_sum += coeff_ak * (v1z + v2z + v3z);
            }
            velocity_field.x() = vx_sum;
            velocity_field.y() = vy_sum;
            velocity_field.z() = vz_sum;

            return velocity_field;
        }

    }
}