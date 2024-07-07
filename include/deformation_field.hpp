#ifndef DEFORMATION_FIELD_HPP
#define DEFORMATION_FIELD_HPP

#include "io.hpp"

#include <eigen3/Eigen/Dense>
#include <cassert>

namespace adi{
    namespace deformation_field{

        struct BasisIndices{
            uint32_t index_1;
            uint32_t index_2;
            uint32_t index_3;
            double eigen_value;
        };
        class DeformationField{
            public:
                DeformationField(const uint32_t &number_of_grid_points) : m_number_of_grid_points(number_of_grid_points){} 

                std::vector<Eigen::Vector3d> computeVelocityField();
                
            private:
                const uint32_t &m_number_of_grid_points;

                /**
                 * @brief Function to compute scalar potential field at every point in the domain
                 * 
                 * @param index_1
                 * @param index_2
                 * @param index3
                 * @param x
                 * @param y
                 * @param z
                 */
                const double phi(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z );
                
                /**
                 * @brief Function to compute partial derivative of scalar potential field with respect to x at every point in the domain
                 * 
                 * @param index_1
                 * @param index_2
                 * @param index3
                 * @param x
                 * @param y
                 * @param z
                 */
                const double dphidx(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z );

                /**
                 * @brief Function to compute partial derivative of scalar potential field with respect to y at every point in the domain
                 * 
                 * @param index_1
                 * @param index_2
                 * @param index3
                 * @param x
                 * @param y
                 * @param z
                 */
                const double dphidy(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z );
                
                /**
                 * @brief Function to compute partial derivative of scalar potential field with respect to z at every point in the domain
                 * 
                 * @param index_1
                 * @param index_2
                 * @param index3
                 * @param x
                 * @param y
                 * @param z
                 */
                const double dphidz(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z );

                const std::vector<BasisIndices> GenerateBasisIndices(const uint32_t &max_number_of_velocity_basis);

                const double ObtainCoefficientOfVelocityField(const double eigen_value);
        };
    }
}
#endif
