#include "deformation_field.hpp"
#include "io.hpp"
#include "constants.hpp"
#include "numerics.hpp"

#include<cmath>


namespace adi{
    namespace deformation_field{
        const double DeformationField::phi(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z )
        {
            return 0.125 * sin(M_PI * index_1 * x) * sin(M_PI * index_2 * y) * sin(M_PI * index_3 * z);
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

        std::tuple<double, double, double> DeformationField::computeVelocityBasisFunctions(const double &coeff_ak, const BasisIndices &base_index, const Eigen::Vector3d &point)
        {
            double vx = 0.0, vy = 0.0, vz = 0.0;

            // Compute partial derivatives of the eigenfunction phi
            double dphi_dx = dphidx(base_index.index_1, base_index.index_2, base_index.index_3, point.x(), point.y(), point.z());
            double dphi_dy = dphidy(base_index.index_1, base_index.index_2, base_index.index_3, point.x(), point.y(), point.z());
            double dphi_dz = dphidz(base_index.index_1, base_index.index_2, base_index.index_3, point.x(), point.y(), point.z());
            

            // Basis function velocities
            double v1x = 0, v1y = dphi_dz, v1z = -dphi_dy;
            double v2x = -dphi_dz, v2y = 0, v2z = dphi_dx;
            double v3x = dphi_dy, v3y = -dphi_dx, v3z = 0;

            vx = coeff_ak * (v1x + v2x + v3x);
            vy = coeff_ak * (v1y + v2y + v3y);
            vz = coeff_ak * (v1z + v2z + v3z);

            return std::make_tuple(vx, vy, vz);
        }

        Eigen::Vector3d DeformationField::computeVelocityField(const Eigen::Vector3d &point, const std::vector<BasisIndices> &base_indices, const Eigen::VectorXd &coeffs_ak)
        {
            Eigen::Vector3d velocity_field;
            std::tuple<double, double, double> velocity_basis_functions;
            double vx_sum = 0.0, vy_sum = 0.0, vz_sum = 0.0;

            for(uint32_t i = 0;i < base_indices.size();++i)
            {
                
                // Compute partial derivatives of the eigenfunction phi
                double coeff_ak;
                if(coeffs_ak.isZero()){
                    coeff_ak = base_indices[i].eigen_value;
                }
                else{
                    coeff_ak = coeffs_ak[i];
                }
                velocity_basis_functions = computeVelocityBasisFunctions(coeff_ak, base_indices.at(i), point);
                vx_sum += std::get<0>(velocity_basis_functions);
                vy_sum += std::get<1>(velocity_basis_functions);
                vz_sum += std::get<2>(velocity_basis_functions);
                
            }
            velocity_field.x() = vx_sum;
            velocity_field.y() = vy_sum;
            velocity_field.z() = vz_sum;

            return velocity_field;
        }

    }
}