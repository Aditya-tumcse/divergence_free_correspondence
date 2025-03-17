#ifndef PARALLEL_OPTIMIZER_HPP
#define PARALLEL_OPTIMIZER_HPP

#include <ceres/ceres.h>
#include <pthread.h>

#include "deformation_field.hpp"
#include "io.hpp"
#include "matching.hpp"
#include "numerics.hpp"
#include "utilities.hpp"

namespace adi {
namespace parallel_optimizer {
class ParallelOptimizer {
public:
  ParallelOptimizer(const unsigned int num_threads)
      : m_num_threads(num_threads) {}

  const unsigned int GetThreadCount() const { return m_num_threads; }

  void Optimize(const std::vector<adi::Point> &target_cloud,
                const std::vector<adi::Point> &src_cloud,
                const Eigen::SparseMatrix<double> &soft_corr_matrix,
                const Eigen::MatrixXd &L_inv,
                const std::array<adi::deformation_field::BasisIndices,
                                 MAX_NUMBER_OF_VELOCITY_BASIS> &base_indices);

  ~ParallelOptimizer() = default;

private:
  const unsigned int m_num_threads;

  struct ThreadWorkSpace {
    size_t start_idx;
    size_t end_idx;
    ceres::Problem *problem;
    const std::vector<adi::Point> *src_cloud;
    const std::vector<adi::Point> *target_cloud;
    const Eigen::SparseMatrix<double> *soft_corr_matrix;
    adi::numerics::PositionIncrementor<double> *position_incrementor;
    const std::array<adi::deformation_field::BasisIndices,
                     MAX_NUMBER_OF_VELOCITY_BASIS> *base_indices;
    pthread_mutex_t *mutex;
  };

  static void *ProcessChunk(void *arg);
};
} // namespace parallel_optimizer
} // namespace adi
#endif