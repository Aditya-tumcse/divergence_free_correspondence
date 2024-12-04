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

const Eigen::MatrixXd computePrecisionMatrix(
    const std::vector<adi::deformation_field::BasisIndices> &basis_indices);

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

const double EucledianDistance(const Eigen::Vector3d &point_1,
                               const Eigen::Vector3d &point_2);

template <typename T> class PositionIncrementor {
public:
  explicit PositionIncrementor(T *const coeffs) : m_coeffs(coeffs) {}

  T *getCoeffs() const { return m_coeffs; }

  T dot(const Eigen::VectorXd &v1) const {
    T sum(0);
    for (int i = 0; i < v1.size(); ++i) {
      sum += m_coeffs[i] * T(v1[i]);
    }
    return sum;
  }

  void computeVelocityField(const Eigen::VectorXd &vel_basis_x,
                            const Eigen::VectorXd &vel_basis_y,
                            const Eigen::VectorXd &vel_basis_z, T *vel_field);

  void incrementPosition(
      T *to_be_updated_pt, const Eigen::VectorXd &vel_basis_x,
      const Eigen::VectorXd &vel_basis_y, const Eigen::VectorXd &vel_basis_z,
      const std::vector<adi::deformation_field::BasisIndices> &basis_indices,
      const double num_time_steps);

  ~PositionIncrementor() = default;

private:
  T *m_coeffs;
};

void Optimize(
    const std::vector<adi::Point> &target_cloud,
    const std::vector<adi::Point> &src_cloud,
    const Eigen::MatrixXd &soft_corr_matrix, const Eigen::MatrixXd &L_inv,
    const std::vector<adi::deformation_field::BasisIndices> &base_indices);
} // namespace numerics
} // namespace adi
#endif