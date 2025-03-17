#include "parallel_optimizer.hpp"

namespace adi {
namespace parallel_optimizer {
void ParallelOptimizer::Optimize(
    const std::vector<adi::Point> &target_cloud,
    const std::vector<adi::Point> &src_cloud,
    const Eigen::SparseMatrix<double> &soft_corr_matrix,
    const Eigen::MatrixXd &L_inv,
    const std::array<adi::deformation_field::BasisIndices,
                     MAX_NUMBER_OF_VELOCITY_BASIS> &base_indices) {
  // Prepare the ceres problem
  ceres::Problem problem;
  double coeffs_ak[MAX_NUMBER_OF_VELOCITY_BASIS];

  // Create local mutex
  pthread_mutex_t problem_mutex;
  pthread_mutex_init(&problem_mutex, nullptr);

  // Initialize the velocity field coefficients with 0
  utilities::fillVector<double>(
      coeffs_ak, Eigen::VectorXd::Zero(MAX_NUMBER_OF_VELOCITY_BASIS));

  // Turn on to initialize the coefficients of velocity field on the basis of
  // the normal distribution of their corresponding eigen values
  // utilities::fillArray<double>(coeffs_ak, base_indices);
  auto position_incrementor =
      std::make_shared<adi::numerics::PositionIncrementor<double>>(coeffs_ak);

  auto workspaces =
      std::make_shared<std::vector<ThreadWorkSpace>>(this->GetThreadCount());
  std::vector<pthread_t> threads(this->GetThreadCount());

  // Divide source cloud amongst threads
  const size_t chunk_size = src_cloud.size() / m_num_threads;

  // Launch the threads
  for (unsigned int tid = 0; tid < this->GetThreadCount(); ++tid) {
    ThreadWorkSpace &thread_ws = (*workspaces)[tid];

    thread_ws.start_idx = tid * chunk_size;
    thread_ws.end_idx =
        (tid == m_num_threads - 1) ? src_cloud.size() : (tid + 1) * chunk_size;

    thread_ws.problem = &problem;
    thread_ws.src_cloud = &src_cloud;
    thread_ws.target_cloud = &target_cloud;
    thread_ws.soft_corr_matrix = &soft_corr_matrix;
    thread_ws.position_incrementor = position_incrementor.get();
    thread_ws.base_indices = &base_indices;
    thread_ws.mutex = &problem_mutex;

    pthread_create(&threads[tid], nullptr, ProcessChunk,
                   static_cast<void *>(&(*workspaces)[tid]));
  }

  for (uint32_t t = 0; t < this->GetThreadCount(); ++t) {
    pthread_join(threads[t], nullptr);
  }

  // Set solver options
  ceres::Solver::Options options;
  options.max_num_iterations = 100;
  options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;

  // Solve the problem
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  // Print solver summary
  std::cout << summary.FullReport() << std::endl;

  pthread_mutex_destroy(&problem_mutex);

  // for (unsigned int id = 0; id < MAX_NUMBER_OF_VELOCITY_BASIS; ++id) {
  //   std::cout << coeffs_ak[id] << std::endl;
  // }
}

void *ParallelOptimizer::ProcessChunk(void *arg) {
  ThreadWorkSpace *tws = static_cast<ThreadWorkSpace *>(arg);
  for (size_t i = tws->start_idx; i < tws->end_idx; ++i) {
    Eigen::Vector3d src_point = tws->src_cloud->at(i).s_point;
    for (size_t j = 0; j < tws->target_cloud->size(); ++j) {
      Eigen::Vector3d target_point = tws->target_cloud->at(j).s_point;

      // TODO: WEIGHTS CAN BE APPLIED TO THE LOSS
      const double weight = tws->soft_corr_matrix->coeff(i, j);
      if (std::isnan(weight))
        continue;

      pthread_mutex_lock(tws->mutex);
      tws->problem->AddResidualBlock(
          new ceres::AutoDiffCostFunction<adi::numerics::MatchingCostFunctor, 1,
                                          MAX_NUMBER_OF_VELOCITY_BASIS>(
              new adi::numerics::MatchingCostFunctor(
                  src_point, target_point, weight, *tws->base_indices)),
          nullptr, tws->position_incrementor->getCoeffs());
      pthread_mutex_unlock(tws->mutex);
    }
  }
  return nullptr;
}
} // namespace parallel_optimizer
} // namespace adi