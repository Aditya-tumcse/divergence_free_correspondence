#ifndef DEFORMATION_FIELD_HPP
#define DEFORMATION_FIELD_HPP

#include "io.hpp"
#include "constants.hpp"
#include "Fastor/Fastor.h"

#include <eigen3/Eigen/Dense>
#include <cassert>

namespace adi{
    namespace deformation_field{
        
        /**
         * @brief Structure for storing the info regarding the basis indices
         */
        struct BasisIndices{
            uint32_t index_1;
            uint32_t index_2;
            uint32_t index_3;
            double eigen_value;
        };
         
        class DeformationField{
            public:
                DeformationField(){}
                
                /**
                 * @brief Computes the deformation field/velocity vector field for the entire pointcloud
                 * 
                 * @param coeffs_ak
                 * @param vel_basis_functions
                 * 
                 * @return An eigen matrix of the velocity field where each column represents each component of the velocity field for each point
                 */
                Eigen::MatrixXd computeVelocityField(const Eigen::VectorXd &coeffs_ak, const Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS, MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH> &vel_basis_functions);  

                // template <typename T>
                // static Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3> computeVelocityFieldTemplated(
                //     const Eigen::Matrix<T, MAX_NUMBER_OF_VELOCITY_BASIS, 1> &coeffs_ak,
                //     const Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS, MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH> &vel_basis_functions);

                /**
                 * @brief Computes velocity basis functions for the pointcloud
                 * 
                 * @param basis_indices
                 * @param points
                 * 
                 * @return Tensor of velocity basis functions for the entire pointcloud and all basis functions
                 */
                Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS, MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH> computeVelocityBasisFunctions(const std::vector<BasisIndices> &basis_indices, const Eigen::MatrixXd &points);
                
                ~DeformationField() = default;
                
            private:
                /**
                 * @brief Function to compute partial derivative of scalar potential field with respect to x at every point in the domain
                 * 
                 * @param basis_indices
                 * @param pts
                 * 
                 * @return Tensor of 1 dimesion that repreents partial derivative of the eigen function with respect to x
                 */
                Fastor::Tensor<double, 1> dphidx(const BasisIndices basis_indices, const Eigen::MatrixXd &pts);

               /**
                 * @brief Function to compute partial derivative of scalar potential field with respect to x at every point in the domain
                 * 
                 * @param basis_indices
                 * @param pts
                 * 
                 * @return Tensor of 1 dimesion that repreents partial derivative of the eigen function with respect to y
                 */
                Fastor::Tensor<double, 1> dphidy(const BasisIndices basis_indices, const Eigen::MatrixXd &pts);
                
                /**
                 * @brief Function to compute partial derivative of scalar potential field with respect to x at every point in the domain
                 * 
                 * @param basis_indices
                 * @param pts
                 * 
                 * @return Tensor of 1 dimesion that repreents partial derivative of the eigen function with respect to z
                 */
                Fastor::Tensor<double, 1> dphidz(const BasisIndices basis_indices, const Eigen::MatrixXd &pts);
        };
    }
}
#endif
