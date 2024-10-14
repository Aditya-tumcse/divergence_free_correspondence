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
                Eigen::Vector3d computeVelocityField(const Eigen::Vector3d &point, const std::vector<BasisIndices> &base_indices, const Eigen::VectorXd &coeffs_ak);  

                std::tuple<double, double, double> computeVelocityBasisFunctions(const double &coeff_ak, const BasisIndices &base_index, const Eigen::Vector3d &point);
                
            private:

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
                
                const double d2phidzdy(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z );

                const double d2phidydz(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z );

                const double d2phidxdz(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z );

                const double d2phidzdx(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z );

                const double d2phidydx(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z );

                const double d2phidxdy(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z );
        };
    }
}
#endif
