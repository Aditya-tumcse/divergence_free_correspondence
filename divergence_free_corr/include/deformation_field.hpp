#ifndef DEFORMATION_FIELD_HPP
#define DEFORMATION_FIELD_HPP

#include "io.hpp"

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <cassert>

namespace adi{
    namespace deformation_field{

        struct BasisIndices{
            uint32_t index_1;
            uint32_t index_2;
            uint32_t index_3;
            double eigen_value;
        };

        void* VelocityBasisWorker(void *arg);
         
        class DeformationField{
            public:
                Eigen::Vector3d computeVelocityField(const Eigen::Vector3d &point, const std::vector<BasisIndices> &base_indices, const Eigen::VectorXd &coeffs_ak);  

                Eigen::Vector3d computeVelocityBasisFunctions(const double &coeff_ak, const BasisIndices &base_index, const Eigen::Vector3d &point);

                Eigen::Tensor<double,3> computeVelocityBasisFunctions(const std::vector<BasisIndices> &basis_indices, const Eigen::MatrixXd &points);
                
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
                // const double dphidx(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z );
                Eigen::Tensor<double, 1> dphidx(const BasisIndices basis_indices, const Eigen::MatrixXd &pts);

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
                Eigen::Tensor<double, 1> dphidy(const BasisIndices basis_indices, const Eigen::MatrixXd &pts);
                
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
                Eigen::Tensor<double, 1> dphidz(const BasisIndices basis_indices, const Eigen::MatrixXd &pts);
        };

        struct VelThreadData{
            uint32_t start_index;
            uint32_t end_index;
            const std::vector<BasisIndices>* basis_indices;
            const Eigen::Vector3d* point;
            const Eigen::VectorXd* coeffs_ak;
            DeformationField* df;
            Eigen::Vector3d* point_vel_field;
        };
    }
}
#endif
