#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include <ceres/ceres.h>
#include <cmath>
#include <eigen3/Eigen/Dense>

#include "constants.hpp"
#include "deformation_field.hpp"
#include "utilities.hpp"

namespace adi {
namespace numerics {
/**
 * @brief Generate vector of deformation field basis indices
 *
 * @param max_number_of_velocity_basis
 *
 * @return Vector of deformation field basis indices
 */
const std::vector<adi::deformation_field::BasisIndices>
GenerateBasisIndices(const uint32_t &max_number_of_velocity_basis);

/**
 * @brief Perform Runge Kutta integration of 2nd order to obtain new positions
 * of the pointcloud to be deformed
 *
 * @param basis_indices
 * @param src_point_cloud
 * @param coeffs_ak
 * @param vel_basis_functions
 * @param num_time_steps
 *
 * @return Eigen matrix of new positions of the pointcloud to be deformed
 */
Eigen::MatrixXd RungeKutta2Integration(
    const std::vector<adi::deformation_field::BasisIndices> &basis_indices,
    Eigen::MatrixXd src_point_cloud, const Eigen::VectorXd &coeffs_ak,
    const Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
                         MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
        &vel_basis_functions,
    const uint32_t num_time_steps);

// template <typename T>
// Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3> RungeKutta2Integration(
//     const std::vector<adi::deformation_field::BasisIndices>& basis_indices,
//     const Eigen::Matrix<T, NUMBER_OF_SAMPLE_POINTS, 3>& src_point_cloud,
//     const Eigen::Matrix<T, MAX_NUMBER_OF_VELOCITY_BASIS, 1>& coeffs_ak,
//     const Fastor::Tensor<T, NUMBER_OF_SAMPLE_POINTS,
//     MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>& vel_basis_functions, const
//     uint32_t num_time_steps);

const double EucledianDistance(const Eigen::Vector3d &point_1,
                               const Eigen::Vector3d &point_2);

template <typename T> class PositionIncrementor {
public:
  explicit PositionIncrementor(Eigen::Vector3<T> *updated_pt)
      : s_updated_pt(updated_pt), s_df() {}

  T *getUpdatedPositions() const { return s_updated_pt; }

  Eigen::Vector3<T> computeVelocityField(const Eigen::VectorX<T> &vel_basis_x,
                                         const Eigen::VectorX<T> &vel_basis_y,
                                         const Eigen::VectorX<T> &vel_basis_z,
                                         const Eigen::VectorX<T> &coeffs);

  void incrementPosition(
      const Eigen::VectorX<T> &vel_basis_x,
      const Eigen::VectorX<T> &vel_basis_y,
      const Eigen::VectorX<T> &vel_basis_z, const Eigen::VectorX<T> &coeffs,
      const std::vector<adi::deformation_field::BasisIndices> &basis_indices,
      const double num_time_steps);

  ~PositionIncrementor() = default;

private:
  Eigen::Vector3<T> *s_updated_pt;
  adi::deformation_field::DeformationField s_df;
};

void Optimize(
    const Fastor::Tensor<double, NUMBER_OF_SAMPLE_POINTS,
                         MAX_NUMBER_OF_VELOCITY_BASIS, TENSOR_DEPTH>
        &vel_basis_functions,
    Eigen::VectorXd *coeffs_ak, const std::vector<adi::Point> &target_cloud,
    const std::vector<adi::Point> &src_cloud,
    const Eigen::MatrixXd &soft_corr_matrix, const Eigen::MatrixXd &L_inv,
    const std::vector<adi::deformation_field::BasisIndices> &base_indices);
} // namespace numerics
} // namespace adi
#endif