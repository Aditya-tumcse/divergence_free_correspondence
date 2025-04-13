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

const std::array<adi::deformation_field::BasisIndices,
                 MAX_NUMBER_OF_VELOCITY_BASIS>
GetUniqueBasisIndices(
    const std::vector<adi::deformation_field::BasisIndices> &basis_indices);

/**
 * @brief Generate vector of deformation field basis indices
 *
 * @param max_number_of_velocity_basis
 *
 * @return Vector of deformation field basis indices
 */
const std::array<adi::deformation_field::BasisIndices,
                 MAX_NUMBER_OF_VELOCITY_BASIS>
GenerateBasisIndices(const uint32_t &max_number_of_velocity_basis);

const Eigen::MatrixXd computePrecisionMatrix(
    const std::array<adi::deformation_field::BasisIndices,
                     MAX_NUMBER_OF_VELOCITY_BASIS> &basis_indices);

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

template <typename T>
T computeResidual(const T *source_pt, const T *target_pt, const T weight) {
  T diff_squared =
      (target_pt[0] - source_pt[0]) * (target_pt[0] - source_pt[0]) +
      (target_pt[1] - source_pt[1]) * (target_pt[1] - source_pt[1]) +
      (target_pt[2] - source_pt[2]) * (target_pt[2] - source_pt[2]);

  T result = weight * ceres::sqrt(diff_squared);

  return result;
}

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
      const std::array<adi::deformation_field::BasisIndices,
                       MAX_NUMBER_OF_VELOCITY_BASIS> &basis_indices,
      const double num_time_steps);

  ~PositionIncrementor() = default;

private:
  T *m_coeffs;
};

struct MatchingCostFunctor {
  MatchingCostFunctor(
      const Eigen::Vector3d &source_point, const Eigen::Vector3d &target_point,
      const double weight,
      const std::array<adi::deformation_field::BasisIndices,
                       MAX_NUMBER_OF_VELOCITY_BASIS> &base_indices)
      : s_target_point(target_point), s_source_point(source_point),
        s_weight(weight), s_base_indices(base_indices) {

    // TODO: compute velocity basis function here for each point
    adi::deformation_field::DeformationField s_df;
    Eigen::MatrixXd vel_basis_functions =
        s_df.computeVelocityBasisFunctionsPerPoint<double>(base_indices,
                                                           source_point);

    s_vel_basis_function_x = vel_basis_functions.col(0);
    s_vel_basis_function_y = vel_basis_functions.col(1);
    s_vel_basis_function_z = vel_basis_functions.col(2);
  }

  template <typename T> bool checkForNaN(const Eigen::VectorX<T> &vec) const {
    for (int i = 0; i < vec.size(); ++i) {
      if (ceres::isnan(vec[i])) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool operator()(const T *const coeffs_ak, T *residual) const {
    T src_pt[3];
    T target_pt[3];
    utilities::toTemplatedPoint(s_source_point, src_pt);
    utilities::toTemplatedPoint(s_target_point, target_pt);

    PositionIncrementor<T> position_incrementor((T *const)(coeffs_ak));
    position_incrementor.incrementPosition(
        src_pt, s_vel_basis_function_x, s_vel_basis_function_y,
        s_vel_basis_function_z, s_base_indices, NUMBER_OF_TIME_STEPS);

    // if (checkForNaN<T>(updated_pos)) {
    //   std::cout << "Updated position is invalid" << std::endl;
    //   return false;
    // }

    residual[0] = computeResidual(src_pt, target_pt, T(s_weight));
    return true;
  }

  const Eigen::Vector3d s_source_point;
  const Eigen::Vector3d s_target_point;
  Eigen::VectorXd s_vel_basis_function_x;
  Eigen::VectorXd s_vel_basis_function_y;
  Eigen::VectorXd s_vel_basis_function_z;
  const double s_weight;
  const std::array<adi::deformation_field::BasisIndices,
                   MAX_NUMBER_OF_VELOCITY_BASIS> &s_base_indices;
};

void Optimize(const std::vector<adi::Point> &target_cloud,
              const std::vector<adi::Point> &src_cloud,
              const Eigen::MatrixXd &soft_corr_matrix,
              const Eigen::MatrixXd &L_inv,
              const std::array<adi::deformation_field::BasisIndices,
                               MAX_NUMBER_OF_VELOCITY_BASIS> &base_indices);
} // namespace numerics
} // namespace adi
#endif